/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_H

/**
 * SolverMVA over a SolverLN layer.
 *
 * Ports of matlab/src/solvers/MVA/solver_mva.m (exact BCMP MVA),
 * solver_amvald.m + solver_amvald_forward.m (the general approximate MVA), and
 * solver_amva.m + solver_mva_analyzer.m (the method dispatch that chooses
 * between them). The chain aggregation they all sit on is in sn_chain.h.
 *
 * SCOPE. The dispatch reproduces the `default` ladder in full, since it is what
 * decides which algorithm each layer gets and therefore what the numbers are.
 * The algorithms behind it are ported for the disciplines a layered model
 * produces -- INF, PS, FCFS, SIRO, LCFSPR -- with the `default` and `seidmann`
 * multiserver approximations and the `default` high-variance setting. Priority
 * (HOL), discriminatory sharing (DPS), load- and class-dependence, the `suri`
 * and `softmin` multiserver variants, the queue-line and fraction-line
 * estimators and open classes are REFUSED by name. Each is a distinct
 * approximation, and mapping any of them onto a neighbour would return a
 * plausible number that is not the reference's.
 *
 * ARITHMETIC. Everything here is field arithmetic. The AMVA loops stop on a
 * tolerance, so an exact run returns the iterate its stopping rule selected,
 * computed without rounding; see the gate note in pfqn_egflinearizer.h.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "line/api/npfqn/npfqn_sqd.h"
#include "line/api/pfqn/pfqn_ab_amva.h"
#include "line/api/pfqn/pfqn_aql.h"
#include "line/api/pfqn/pfqn_qsa.h"
#include "line/api/pfqn/pfqn_lcfsqn_mva.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_chow.h"
#include "line/api/pfqn/pfqn_clust.h"
#include "line/api/pfqn/pfqn_dmlin.h"
#include "line/api/pfqn/pfqn_lcp.h"
#include "line/api/pfqn/pfqn_pam.h"
#include "line/api/pfqn/pfqn_cdfun.h"
#include "line/api/pfqn/pfqn_conwayms.h"
#include "line/api/pfqn/pfqn_jdfun.h"
#include "line/api/pfqn/pfqn_lldfun.h"
#include "line/api/pfqn/pfqn_scat.h"
#include "line/api/pfqn/pfqn_tay.h"
#include "line/api/pfqn/pfqn_linearizermx.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/api/pfqn/pfqn_schmidt.h"
#include "line/api/pfqn/pfqn_schmidt_ext.h"
#include "line/api/pfqn/pfqn_sqni.h"
#include "line/api/sum/sum_closed.h"
#include "line/api/sum/sum_closing.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"

namespace line {
namespace mva {

// ---------------------------------------------------------------------------
// Exact MVA
// ---------------------------------------------------------------------------

/**
 * Port of `solver_mva_lcfsqn.m`: the closed two-station LCFS + LCFS-PR network.
 *
 * `pfqn_lcfsqn_mva` returns the throughput, the queue lengths and the
 * utilizations for the pair directly; everything else here is the mapping back
 * onto the model's station indexing, with the response times from Little's law
 * and the cycle time as their sum. `lG` is NaN: this recursion carries no
 * normalizing constant, and the reference says so rather than reporting zero.
 */
template <class T>
MvaSolution<T> solver_mva_lcfsqn(const qn::NetworkStruct<T>& L, std::size_t lcfs_ist,
                                 std::size_t lcfspr_ist) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, R = L.nclasses;

    // A CLASS OF ZERO POPULATION NEEDS A PLACEHOLDER RATE, not a zero. The
    // reference initializes these to zeros (solver_mva_lcfsqn.m:31-34) and gets
    // away with it because pfqn_lcfsqn_mva.m validates nothing; the C++ api
    // added a check that rejects a non-positive alpha for EVERY class
    // regardless of its population (pfqn_lcfsqn_mva.h:97-99), so a closed LCFS
    // model carrying an empty class threw where MATLAB answers. One is inert:
    // it enters only as 1^0 in the product A = prod alpha_r^{n_r}, as
    // alpha_r * T_r with T_r = 0 in the throughput sum, and as U = T*alpha = 0.
    std::vector<T> alpha(R, one), beta(R, one);
    std::vector<int> N(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        const double nr = L.classes[r].population;
        if (!(nr > 0.0)) continue;
        N[r] = static_cast<int>(std::llround(nr));
        const T mu_l = L.rates(lcfs_ist - 1, r);
        const T mu_p = L.rates(lcfspr_ist - 1, r);
        if (!(mu_l > zero) || !std::isfinite(num_traits<T>::to_double(mu_l)))
            throw InputError("solver_mva_lcfsqn: invalid service rate at the LCFS station for class " +
                             std::to_string(r + 1));
        if (!(mu_p > zero) || !std::isfinite(num_traits<T>::to_double(mu_p)))
            throw InputError("solver_mva_lcfsqn: invalid service rate at the LCFS-PR station for "
                             "class " + std::to_string(r + 1));
        alpha[r] = T(one / mu_l);
        beta[r] = T(one / mu_p);
    }

    const pfqn::LcfsMvaResult<T> res = pfqn::pfqn_lcfsqn_mva(alpha, beta, N);

    MvaSolution<T> out;
    out.Q = Matrix<T>(M, R, zero);
    out.U = Matrix<T>(M, R, zero);
    out.R = Matrix<T>(M, R, zero);
    out.Tp = Matrix<T>(M, R, zero);
    out.C.assign(R, zero);
    out.X.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        out.Q(lcfs_ist - 1, r) = res.Q(0, r);
        out.Q(lcfspr_ist - 1, r) = res.Q(1, r);
        out.U(lcfs_ist - 1, r) = res.U(0, r);
        out.U(lcfspr_ist - 1, r) = res.U(1, r);
        if (N[r] <= 0) continue;
        out.X[r] = res.T_[r];
        out.Tp(lcfs_ist - 1, r) = res.T_[r];
        out.Tp(lcfspr_ist - 1, r) = res.T_[r];
    }
    for (std::size_t i : {lcfs_ist - 1, lcfspr_ist - 1})
        for (std::size_t r = 0; r < R; ++r)
            if (out.Tp(i, r) > zero) out.R(i, r) = T(out.Q(i, r) / out.Tp(i, r));
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > 0) out.C[r] = T(out.R(lcfs_ist - 1, r) + out.R(lcfspr_ist - 1, r));
    // REPORTED AS 'exact', not as 'lcfsqn'. `solver_mva_analyzer` names the
    // method AFTER solver_mva returns, so the reference reports every internal
    // branch of it under the one name; a caller that saw 'lcfsqn' here would
    // disagree with MATLAB on `actualmethod` while agreeing on every number.
    out.method = "exact";
    out.lG = std::numeric_limits<double>::quiet_NaN();
    return out;
}

/** Port of solver_mva.m, closed and mixed product-form networks. */
template <class T>
MvaSolution<T> solver_mva(const qn::NetworkStruct<T>& L, const ChainDemands<T>& d, const MvaOptions& opt) {
    using qn::SchedStrategy;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, C = L.nchains;

    // The special-cased LCFS + LCFS-PR pair, BEFORE the product-form test: an
    // LCFS station is not product-form on its own, and the reference reaches
    // this branch on a model the generic path would reject.
    std::vector<std::size_t> lcfs, lcfspr;  // 1-based station indices
    for (std::size_t i = 0; i < M; ++i) {
        if (L.stations[i].sched == SchedStrategy::LCFS) lcfs.push_back(i + 1);
        if (L.stations[i].sched == SchedStrategy::LCFSPR) lcfspr.push_back(i + 1);
    }
    if (!lcfs.empty() && !lcfspr.empty()) {
        if (lcfs.size() != 1 || lcfspr.size() != 1)
            throw UnsupportedError(
                "solver_mva: LCFS MVA requires exactly one LCFS and one LCFS-PR station");
        for (std::size_t c = 0; c < C; ++c)
            if (std::isinf(d.Nchain[c]))
                throw UnsupportedError("solver_mva: LCFS MVA requires a closed queueing network");
        // A self-loop would let a job re-enter the station it just left, which
        // the two-station recursion has no term for.
        const std::size_t S = L.nof_stateful(), Rn = L.nclasses;
        for (std::size_t ist : {lcfs[0], lcfspr[0]}) {
            const std::size_t sf = L.stateful_of_station(ist) - 1;
            for (std::size_t r = 0; r < Rn; ++r)
                if (L.rt.rows() == S * Rn && L.rt(sf * Rn + r, sf * Rn + r) > zero)
                    throw UnsupportedError(
                        "solver_mva: LCFS MVA does not support self-loops at stations");
        }
        return solver_mva_lcfsqn(L, lcfs[0], lcfspr[0]);
    }
    if (!lcfs.empty())
        throw UnsupportedError("solver_mva: LCFS scheduling requires a paired LCFS-PR station");

    // METHOD 'mva' IS THE DELIBERATE APPROXIMATION: the dispatch warns that the
    // exact recursion is being run outside its hypotheses and promises an answer,
    // so throwing here would contradict its own message. Only an implicit or
    // 'exact' request is refused.
    if (!L.has_product_form() && opt.method != "mva")
        throw UnsupportedError("solver_mva: the layer does not have a product form");

    std::vector<std::size_t> infSET, qSET;  // 0-based station indices
    for (std::size_t i = 0; i < M; ++i) {
        switch (L.stations[i].sched) {
            case SchedStrategy::EXT: break;
            case SchedStrategy::INF: infSET.push_back(i); break;
            case SchedStrategy::PS:
            case SchedStrategy::LCFSPR:
            case SchedStrategy::FCFS:
            case SchedStrategy::SIRO: qSET.push_back(i); break;
            default:
                throw UnsupportedError(std::string("solver_mva: unsupported exact MVA analysis for ") +
                                       lang::sched_to_text(L.stations[i].sched) + " scheduling");
        }
    }

    // demands at the queueing and delay stations, chain by chain
    Matrix<T> Lq(qSET.size(), C, zero), Zd(infSET.size(), C, zero);
    for (std::size_t a = 0; a < qSET.size(); ++a)
        for (std::size_t c = 0; c < C; ++c)
            Lq(a, c) = T(d.STchain(qSET[a], c) * d.Vchain(qSET[a], c));
    for (std::size_t a = 0; a < infSET.size(); ++a)
        for (std::size_t c = 0; c < C; ++c)
            Zd(a, c) = T(d.STchain(infSET[a], c) * d.Vchain(infSET[a], c));

    // open-chain arrival-rate rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    std::vector<T> lambda(C, zero);
    std::vector<int> N(C, 0);
    for (std::size_t c = 0; c < C; ++c) {
        if (!std::isfinite(d.Nchain[c])) {
            N[c] = pfqn::OPEN_CLASS;
            const T st = d.STchain(d.refstatchain[c] - 1, c);
            lambda[c] = st > zero ? T(num_traits<T>::from_int(1) / st) : zero;
            continue;
        }
        N[c] = static_cast<int>(std::llround(d.Nchain[c]));
    }
    std::vector<int> S;
    for (std::size_t a = 0; a < qSET.size(); ++a) {
        const double s = L.stations[qSET[a]].nservers;
        if (!std::isfinite(s))
            throw UnsupportedError("solver_mva: a queueing station has infinitely many servers but "
                                   "is not inf-scheduled");
        S.push_back(static_cast<int>(std::llround(s)));
    }

    // Interlocked flow (Franks 1999, Eq. 4.7): a request cannot queue behind work that
    // its own submission caused, so the arrival-instant queue drops the interlocked share
    // of the other chains. SolverLN supplies the matrix, class-indexed.
    const Matrix<T> IL = sn_interlock_chain<T>(L, opt.interlock);
    // the interlocked recursion is a separate entry point: pfqn_mvams and the pfqn_mva
    // family it dispatches to carry the standard arrival theorem only
    const pfqn::MvaResult<T> pf =
        IL.empty() ? pfqn::pfqn_mvams(lambda, Lq, N, Zd, std::vector<int>(qSET.size(), 1), S)
                   : pfqn::pfqn_mvams_ilock(lambda, Lq, N, Zd, std::vector<int>(qSET.size(), 1), S, IL);

    Matrix<T> Qchain(M, C, zero), Wchain(M, C, zero), Tchain(M, C, zero), Uchain(M, C, zero);
    std::vector<T> Xchain = pf.XN;
    for (std::size_t a = 0; a < qSET.size(); ++a)
        for (std::size_t c = 0; c < C; ++c) Qchain(qSET[a], c) = pf.QN(a, c);
    for (std::size_t a = 0; a < infSET.size(); ++a)
        for (std::size_t c = 0; c < C; ++c)
            Qchain(infSET[a], c) = T(Xchain[c] * d.STchain(infSET[a], c) * d.Vchain(infSET[a], c));

    std::vector<std::size_t> rset;
    for (std::size_t c = 0; c < C; ++c)
        if (d.Nchain[c] != 0.0) rset.push_back(c);

    for (std::size_t c : rset) {
        for (std::size_t i : infSET) Wchain(i, c) = d.STchain(i, c);
        for (std::size_t i : qSET) {
            if (std::isinf(L.stations[i].nservers)) {
                Wchain(i, c) = d.STchain(i, c);
            } else if (d.Vchain(i, c) == zero || Xchain[c] == zero) {
                Wchain(i, c) = zero;
            } else {
                Wchain(i, c) = T(Qchain(i, c) / (Xchain[c] * d.Vchain(i, c)));
            }
        }
    }

    std::vector<T> Cc(C, zero);
    for (std::size_t c : rset) {
        T sw = zero;
        for (std::size_t i = 0; i < M; ++i) sw += Wchain(i, c);
        if (sw == zero) {
            Xchain[c] = zero;
        } else {
            T cyc = zero;
            for (std::size_t i = 0; i < M; ++i) cyc += d.Vchain(i, c) * Wchain(i, c);
            Cc[c] = cyc;
            // open-chain throughput rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
            if (cyc != zero && std::isfinite(d.Nchain[c]))
                Xchain[c] = T(num_traits<T>::from_double(d.Nchain[c]) / cyc);
        }
        for (std::size_t i = 0; i < M; ++i) {
            Qchain(i, c) = T(Xchain[c] * d.Vchain(i, c) * Wchain(i, c));
            Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c : rset) {
            const T u = T(d.Vchain(i, c) * d.STchain(i, c) * Xchain[c]);
            Uchain(i, c) = std::isinf(L.stations[i].nservers)
                               ? u
                               : T(u / num_traits<T>::from_double(L.stations[i].nservers));
        }

    // renormalise a utilization that the product-form formula pushed past one
    for (std::size_t i = 0; i < M; ++i) {
        if (!(L.stations[i].sched == SchedStrategy::FCFS || L.stations[i].sched == SchedStrategy::PS))
            continue;
        T usum = zero;
        for (std::size_t c = 0; c < C; ++c) usum += Uchain(i, c);
        if (num_traits<T>::to_double(usum) <= 1.0 + opt.tol) continue;
        T den = zero;
        for (std::size_t c = 0; c < C; ++c) den += d.Vchain(i, c) * d.STchain(i, c) * Xchain[c];
        if (den == zero) continue;
        const T cap = usum > num_traits<T>::from_int(1) ? num_traits<T>::from_int(1) : usum;
        for (std::size_t c = 0; c < C; ++c) {
            if (!(num_traits<T>::to_double(d.Vchain(i, c) * d.STchain(i, c)) > opt.tol)) continue;
            Uchain(i, c) = T(cap * d.Vchain(i, c) * d.STchain(i, c) * Xchain[c] / den);
        }
    }

    Matrix<T> Rchain(M, C, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c)
            if (Tchain(i, c) != zero) Rchain(i, c) = T(Qchain(i, c) / Tchain(i, c));
    for (std::size_t c = 0; c < C; ++c) {
        if (d.Nchain[c] != 0.0) continue;
        Xchain[c] = zero;
        for (std::size_t i = 0; i < M; ++i) {
            Uchain(i, c) = zero;
            Qchain(i, c) = zero;
            Rchain(i, c) = zero;
            Tchain(i, c) = zero;
            Wchain(i, c) = zero;
        }
    }

    const ClassResults<T> cr = sn_deaggregate_chain_results(
        L, d, Matrix<T>(), Matrix<T>(), Rchain, Tchain, Xchain);
    MvaSolution<T> out;
    out.Q = cr.Q;
    out.U = cr.U;
    out.R = cr.R;
    out.Tp = cr.Tp;
    out.C = cr.C;
    out.X = cr.X;
    out.method = "exact";
    out.lG = pf.lG;
    return out;
}

// ---------------------------------------------------------------------------
// Approximate MVA, the general (chain-level) path
// ---------------------------------------------------------------------------

namespace detail {

/**
 * The multiserver capacity term of solver_amvald_forward, all four rules.
 *
 *   softmin    1/softmin(n_i, c_i, 20) at EVERY station
 *   seidmann   1/c_i, and 1 at an infinite server
 *   default    softmin, then the FCFS-family stations overridden with 1/c_i
 *   suri       1 everywhere; the multiplicity is carried by the Suri factor
 *              applied to the queue length inside the residence formula
 *
 * softmin is a genuine transcendental (exp), so the exact backend can evaluate
 * this only where softmin is not reached: seidmann and suri everywhere, and
 * default at an infinite server or an FCFS-family station. Anything else
 * refuses by name.
 */
template <class T>
std::vector<T> ms_term(const qn::NetworkStruct<T>& L, const std::vector<T>& narrival,
                       const std::string& multiserver) {
    using qn::SchedStrategy;
    const std::size_t M = L.nstations;
    const T one = num_traits<T>::from_int(1);
    std::vector<T> r(M, one);
    if (multiserver == "suri") return r;  // no demand scaling; see suri_factor
    if (!(multiserver == "default" || multiserver == "seidmann" || multiserver == "softmin"))
        throw UnsupportedError("solver_amvald: multiserver approximation '" + multiserver +
                               "' is not implemented in this port");
    for (std::size_t i = 0; i < M; ++i) {
        const double c = L.stations[i].nservers;
        const SchedStrategy sc = L.stations[i].sched;
        const bool fcfs_family =
            sc == SchedStrategy::FCFS || sc == SchedStrategy::SIRO || sc == SchedStrategy::LCFSPR;
        if (multiserver == "seidmann" || (multiserver == "default" && fcfs_family)) {
            r[i] = std::isinf(c) ? one : T(one / num_traits<T>::from_double(c));
            continue;
        }
        if (std::isinf(c)) {
            r[i] = one;  // delay server, pfqn_lldfun returns 1
            continue;
        }
        if constexpr (num_traits<T>::has_transcendental) {
            // 1/softmin(n, c, alpha), in the reference's numerically safe form
            const double alpha = 20.0;
            const double n = num_traits<T>::to_double(narrival[i]);
            const double lo = std::min(n, c), hi = std::max(n, c);
            const double gap = hi - lo;
            const double w = alpha * gap > 745.0 ? 0.0 : std::exp(-alpha * gap);
            const double sm = lo + gap * w / (1.0 + w);
            r[i] = num_traits<T>::from_double(sm > 0.0 ? 1.0 / sm : 1.0 / std::min(n, c));
        } else {
            throw UnsupportedError(
                "solver_amvald: the soft-minimum multiserver term evaluates exp, which exact "
                "arithmetic has no representation for; use the double or real backend, the "
                "seidmann or suri multiserver rule, or a layer whose queueing stations are "
                "FCFS-family or infinite-server");
        }
    }
    return r;
}

/**
 * The Suri multiserver factor, applied to the arrival queue length rather than
 * to the demand: rho^(4.464 (c^0.676 - 1)) / c at a finite multiserver station,
 * 1/c when its utilization is zero, and 0 at an infinite or zero server.
 *
 * The utilization it reads is the ITERATE's, so the factor moves with the fixed
 * point; the constants are the reference's fitted ones.
 */
template <class T>
std::vector<T> suri_factor(const qn::NetworkStruct<T>& L, const Matrix<T>& Uin, double tol) {
    const std::size_t M = L.nstations, K = Uin.cols();
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    std::vector<T> f(M, one);
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_amvald: the Suri multiserver factor is a real power of the utilization, "
            "which exact arithmetic has no representation for; use the double or real backend");
    } else {
        const double alpha = 4.464, beta = 0.676;
        for (std::size_t i = 0; i < M; ++i) {
            const double c = L.stations[i].nservers;
            if (std::isinf(c) || c == 0.0) {
                f[i] = zero;
                continue;
            }
            if (!(c > 1.0)) continue;  // single server: the factor stays 1
            T usum = zero;
            for (std::size_t s = 0; s < K; ++s) usum += Uin(i, s);
            double rho = num_traits<T>::to_double(usum) / c;
            if (rho > 1.0 - tol) rho = 1.0 - tol;
            f[i] = rho > 0.0 ? num_traits<T>::from_double(
                                   std::pow(rho, alpha * (std::pow(c, beta) - 1.0)) / c)
                             : T(one / num_traits<T>::from_double(c));
        }
        return f;
    }
}

}  // namespace detail

/**
 * Port of solver_amvald.m together with solver_amvald_forward.m, restricted to
 * the INF / PS / FCFS-family disciplines with no load or class dependence.
 */
template <class T>
MvaSolution<T> solver_amvald(const qn::NetworkStruct<T>& L, const ChainDemands<T>& d,
                             const MvaOptions& opt, const std::string& method, bool& converged,
                             const Matrix<T>& init_sol = Matrix<T>()) {
    using qn::SchedStrategy;
    // EXACT ARITHMETIC HAS NO BOUNDED ITERATE HERE, which is why this refuses
    // rather than runs slowly. The exact product-form recursion is a fixed
    // number of field operations and stays rational; this is a Picard fixed
    // point, and every sweep multiplies the denominators of the one before it,
    // so the bit length grows geometrically in the iteration count. A layer of
    // lqn_ofbiz that solves in milliseconds under `double` was measured at over
    // four and a half hours under `Rational`, with operands millions of bits
    // wide and no answer -- and it only reaches here at all because the model
    // failed the product-form gate. Refusing by name is the same rule the
    // load-dependent soft minimum and the Suri factor follow below.
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_amvald: the approximate MVA is an iterative fixed point whose iterates "
            "accumulate the product of every denominator seen so far, so exact rational "
            "arithmetic grows without bound and the solve does not terminate; rerun with "
            "--arith double or --arith real, or give the model a product form (BCMP type 1 "
            "asks FCFS service to be exponential) so the exact recursion applies");
    }
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nchains;

    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy s = L.stations[i].sched;
        if (!(s == SchedStrategy::INF || s == SchedStrategy::PS || s == SchedStrategy::FCFS ||
              s == SchedStrategy::SIRO || s == SchedStrategy::LCFSPR || s == SchedStrategy::EXT ||
              s == SchedStrategy::DPS || s == SchedStrategy::HOL ||
              s == SchedStrategy::FCFSPRPRIO))
            throw UnsupportedError(std::string("solver_amvald: ") + lang::sched_to_text(s) +
                                   " scheduling is not implemented in this port");
    }
    // A PURE-DELAY MODEL SELECTS NO APPROXIMATION, so the name is not asked
    // about there. With no queueing station the arrival-instant queue length is
    // identically zero, every AMVA scheme collapses to the same exact delay
    // solution, and this analyzer is where `solver_amva` sends that model
    // BECAUSE the product-form kernels would be handed a zero-row demand matrix
    // (solver_amva.m says so in as many words). Asking the whitelist first
    // refused twelve names -- ab, schmidt, schmidt-ext, scat, lcp, chow, pamb,
    // pami, pamt, clust, dmlin and qsa -- on a model that needs no arm for any
    // of them, while the report went on offering all twelve.
    bool any_queue = false;
    for (std::size_t i = 0; i < M && !any_queue; ++i)
        if (L.stations[i].sched != SchedStrategy::INF && L.stations[i].sched != SchedStrategy::EXT)
            any_queue = true;
    const bool linmethod = (method == "lin" || method == "qdlin");
    if (any_queue &&
        !(linmethod || method == "qd" || method == "default" || method == "bs" ||
          method == "egflin" || method == "gflin" || method == "qli" || method == "fli" ||
          method == "aql" || method == "qdaql" || method == "tay" || method == "priomva"))
        throw UnsupportedError("solver_amvald: method '" + method + "' is not implemented");

    // Nt closed-population-only rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    double Nt = 0.0;
    for (std::size_t c = 0; c < K; ++c)
        if (std::isfinite(d.Nchain[c])) Nt += d.Nchain[c];
    const T deltaT = Nt > 0.0 ? num_traits<T>::from_double((Nt - 1.0) / Nt) : one;
    // deltaclass open-chain rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    std::vector<T> deltaclass(K, one);
    for (std::size_t c = 0; c < K; ++c)
        if (std::isfinite(d.Nchain[c]) && d.Nchain[c] > 0.0)
            deltaclass[c] = num_traits<T>::from_double((d.Nchain[c] - 1.0) / d.Nchain[c]);

    // nnz/ccl/ocl rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    std::vector<std::size_t> nnz, ccl, ocl;
    std::vector<bool> isopen(K, false);
    for (std::size_t c = 0; c < K; ++c) {
        if (!(d.Nchain[c] > 0.0)) continue;
        nnz.push_back(c);
        if (std::isfinite(d.Nchain[c])) ccl.push_back(c);
        else { ocl.push_back(c); isopen[c] = true; }
    }

    // warm-start-under-outer-fixed-point rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    Matrix<T> Qchain(M, K, zero);
    if (init_sol.rows() == M && init_sol.cols() == K) {
        Qchain = init_sol;
    } else {
        // open-chain warm-start zeroing rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
        for (std::size_t c = 0; c < K; ++c) {
            if (!std::isfinite(d.Nchain[c])) continue;
            for (std::size_t i = 0; i < M; ++i)
                Qchain(i, c) =
                    T(num_traits<T>::from_double(d.Nchain[c]) / num_traits<T>::from_int(int(M)));
        }
    }

    std::vector<T> Xchain(K, zero);
    for (std::size_t c = 0; c < K; ++c) {
        if (isopen[c]) {
            // open-chain arrival-rate rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
            const T st = d.STchain(d.refstatchain[c] - 1, c);
            Xchain[c] = st > zero ? T(one / st) : zero;
            continue;
        }
        T s = zero;
        for (std::size_t i = 0; i < M; ++i) s += d.STchain(i, c);
        if (s != zero) Xchain[c] = T(one / s);
    }
    Matrix<T> Uchain(M, K, zero), Tchain(M, K, zero), Wchain(M, K, zero), STeff = d.STchain;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c : nnz) {
            const T u = T(d.Vchain(i, c) * d.STchain(i, c) * Xchain[c]);
            Uchain(i, c) = std::isinf(L.stations[i].nservers)
                               ? u
                               : T(u / num_traits<T>::from_double(L.stations[i].nservers));
        }

    // gamma(s, i, r) for 'lin'; gamma(s, i) otherwise
    std::vector<Matrix<T>> gamma(K, Matrix<T>(M, K, zero));

    const T omicron = num_traits<T>::from_double(0.5);
    int totiter = 0;
    const int max_totiter = std::min(opt.iter_max, 10000);
    line::util::LineConsole::loop("running the AMVA fixed point (tolerance %g)", opt.tol);
    const int inner_cap = static_cast<int>(std::sqrt(static_cast<double>(opt.iter_max)));

    // priority set rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    const bool has_prio = [&] {
        int lo = 0, hi = 0;
        bool first = true;
        for (std::size_t c = 0; c < K; ++c) {
            const int p = L.classes.empty() ? 0 : L.classes[c].prio;
            if (first) { lo = hi = p; first = false; }
            lo = std::min(lo, p);
            hi = std::max(hi, p);
        }
        return hi != lo;
    }();
    std::vector<std::vector<std::size_t>> ehprio(K), hprio(K), eprio(K), lprio(K);
    for (std::size_t r : nnz) {
        if (!has_prio) {
            // every class ties, so HOL degenerates to FCFS: all classes are equal-or-higher
            for (std::size_t s : nnz) {
                ehprio[r].push_back(s);
                eprio[r].push_back(s);
            }
            continue;
        }
        for (std::size_t s : nnz) {
            const int pr = L.classes[r].prio, ps = L.classes[s].prio;
            if (ps <= pr) ehprio[r].push_back(s);
            if (ps < pr) hprio[r].push_back(s);
            if (ps == pr) eprio[r].push_back(s);
            if (ps > pr) lprio[r].push_back(s);
        }
    }

    // lldscaling/cdscaling emptiness rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    std::size_t smax = 0;
    for (const auto& st : L.stations) smax = std::max(smax, st.lldscaling.size());
    Matrix<T> lldscaling;
    if (smax > 0) {
        lldscaling = Matrix<T>(M, smax, one);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < L.stations[i].lldscaling.size(); ++k)
                lldscaling(i, k) = L.stations[i].lldscaling[k];
    }
    std::vector<lang::CdScaling<T>> cdscaling;
    for (const auto& st : L.stations)
        if (st.cdscaling) {
            cdscaling.assign(M, lang::CdScaling<T>());
            break;
        }
    if (!cdscaling.empty())
        for (std::size_t i = 0; i < M; ++i) cdscaling[i] = L.stations[i].cdscaling;
    // sn.jdscaling is a SEPARATE list, never folded into cdscaling: the two are
    // evaluated at DIFFERENT points below, so a single list would lose the one
    // piece of information that tells them apart. Assembled by the same rule --
    // empty means "no station has one", which is what pfqn_jdfun tests.
    std::vector<lang::CdScaling<T>> jdscaling;
    for (const auto& st : L.stations)
        if (st.jdscaling) {
            jdscaling.assign(M, lang::CdScaling<T>());
            break;
        }
    if (!jdscaling.empty())
        for (std::size_t i = 0; i < M; ++i) jdscaling[i] = L.stations[i].jdscaling;
    const bool has_scaling = (smax > 0) || !cdscaling.empty() || !jdscaling.empty();

    // tau(s,r) rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    std::vector<std::vector<T>> tau(K, std::vector<T>(K, zero));
    // Throughput vector the tau differences were taken against, so that Xref[r] + tau[s][r]
    // is an arrival-instant throughput from one and the same sweep. Adding tau to the moving
    // inner iterate instead mixes two sweeps and can exceed the service capacity. Empty when
    // no Linearizer recursion runs, and then the inner iterate is the reference.
    std::vector<T> Xref;

    // Interlocked flow (Franks 1999, Eq. 4.7): the option carries the matrix CLASS-indexed,
    // and the forward step below works in the chain basis.
    const Matrix<T> ILchain = sn_interlock_chain<T>(L, opt.interlock);

    // one forward evaluation: residence times from the current queue lengths
    auto forward = [&](const Matrix<T>& Qin, const std::vector<T>& Xin, const Matrix<T>& Uin,
                       const std::vector<double>& Nc, Matrix<T>& Wout, Matrix<T>& STeffOut) {
        double Ntl = 0.0;
        for (std::size_t c = 0; c < K; ++c)
            if (std::isfinite(Nc[c])) Ntl += Nc[c];
        const T deltaL = Ntl > 0.0 ? num_traits<T>::from_double((Ntl - 1.0) / Ntl) : one;
        std::vector<T> dcl(K, one);
        for (std::size_t c = 0; c < K; ++c)
            if (std::isfinite(Nc[c]) && Nc[c] > 0.0)
                dcl[c] = num_traits<T>::from_double((Nc[c] - 1.0) / Nc[c]);

        std::vector<T> interpTot(M, zero);
        Matrix<T> selfArvl(M, K, zero), totArvl(M, K, zero), totArvlOpen(M, K, zero);
        for (std::size_t i = 0; i < M; ++i) {
            T s = zero;
            for (std::size_t c : nnz) s += Qin(i, c);
            interpTot[i] = T(deltaL * s);
            const bool hol = L.stations[i].sched == SchedStrategy::HOL;
            for (std::size_t c : nnz) {
                selfArvl(i, c) = T(dcl[c] * Qin(i, c));
                if (hol) {
                    T eh = zero, ehx = zero;
                    for (std::size_t s2 : ehprio[c]) {
                        eh += Qin(i, s2);
                        if (s2 != c) ehx += Qin(i, s2);
                    }
                    totArvlOpen(i, c) = eh;
                    totArvl(i, c) = T(dcl[c] * Qin(i, c) + ehx);
                } else {
                    totArvlOpen(i, c) = s;
                    totArvl(i, c) = T(dcl[c] * Qin(i, c) + s - Qin(i, c));
                }
            }
        }

        // multiserver softmin population rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
        std::vector<T> gmean(M, zero);
        if (!ccl.empty()) {
            if (linmethod) {
                // g(s,i) = sum_r deltaL * N_r * gamma(s,i,r); then mean over s
                const T nccl = num_traits<T>::from_int(int(ccl.size()));
                for (std::size_t i = 0; i < M; ++i) {
                    T acc = zero;
                    for (std::size_t s : ccl) {
                        T gs = zero;
                        for (std::size_t r : ccl)
                            gs += T(num_traits<T>::from_double(Nc[r]) * gamma[s](i, r));
                        acc += T(deltaL * gs);
                    }
                    gmean[i] = T(acc / nccl);
                }
            } else {
                // mean(g) scalar-collapse rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                T acc = zero;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r : ccl)
                        acc += T(num_traits<T>::from_double(Ntl - 1.0) * gamma[r](i, r));
                const T sc = T(acc / num_traits<T>::from_int(int(M)));
                for (std::size_t i = 0; i < M; ++i) gmean[i] = sc;
            }
        }
        std::vector<T> narrival(M, zero);
        for (std::size_t i = 0; i < M; ++i) narrival[i] = T(one + interpTot[i] + gmean[i]);
        const std::vector<T> msterm = detail::ms_term(L, narrival, opt.multiserver);
        std::vector<T> suri(M, one);
        if (opt.multiserver == "suri") suri = detail::suri_factor(L, Uin, opt.tol);

        // The load-dependent term, per class for the Linearizer (whose
        // correction is class specific) and shared otherwise.
        Matrix<T> lldterm(M, K, one);
        if (smax > 0) {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "solver_amvald: load-dependent scaling interpolates a rate lattice with a "
                    "soft minimum, which exact arithmetic has no representation for; use the "
                    "double or real backend");
            } else {
                if (linmethod && !ccl.empty() && !nnz.empty()) {
                    for (std::size_t r : nnz) {
                        std::vector<T> arg(M, zero);
                        for (std::size_t i = 0; i < M; ++i) {
                            T gcorr = zero;
                            for (std::size_t s : ccl)
                                gcorr += T(num_traits<T>::from_double(Nc[s]) * gamma[r](i, s));
                            gcorr -= gamma[r](i, r);
                            arg[i] = T(one + interpTot[i] + gcorr);
                        }
                        const std::vector<T> v =
                            pfqn::pfqn_lldfun(arg, lldscaling, std::vector<double>());
                        for (std::size_t i = 0; i < M; ++i) lldterm(i, r) = v[i];
                    }
                } else {
                    std::vector<T> arg(M, zero);
                    for (std::size_t i = 0; i < M; ++i) arg[i] = T(one + interpTot[i]);
                    const std::vector<T> v =
                        pfqn::pfqn_lldfun(arg, lldscaling, std::vector<double>());
                    for (std::size_t i = 0; i < M; ++i)
                        for (std::size_t c = 0; c < K; ++c) lldterm(i, c) = v[i];
                }
            }
        }

        // class-dependent term rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
        Matrix<T> cdterm(M, K, one);
        if (!cdscaling.empty()) {
            for (std::size_t r : nnz) {
                Matrix<T> arg(M, K, zero);
                const bool closed = std::isfinite(Nc[r]);
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t c = 0; c < K; ++c)
                        arg(i, c) = T(one + (closed ? selfArvl(i, c) : Qin(i, c)));
                if (closed && linmethod)
                    for (std::size_t i = 0; i < M; ++i) {
                        const T gself =
                            T(num_traits<T>::from_double(Nc[r] - 1.0) * gamma[r](i, r));
                        for (std::size_t c = 0; c < K; ++c) arg(i, c) = T(arg(i, c) + gself);
                    }
                const std::vector<T> v = pfqn::pfqn_cdfun(arg, cdscaling, r);
                for (std::size_t i = 0; i < M; ++i) cdterm(i, r) = v[i];
            }
        }

        // joint-dependence term eta_i (non-product-form). It is NOT evaluated at
        // the cdterm point: beta_{i,r} reads only its OWN marginal, so the +1
        // cdterm puts on the non-arriving classes is inert there, while eta reads
        // the WHOLE occupancy row and every coordinate matters. The arrival
        // theorem gives the arriving class-r job one extra job OF ITS OWN CLASS
        // and leaves the others at their means, so the point is
        //   eta_i(Q_{i,1}, ..., 1 + delta_r Q_{i,r}, ..., Q_{i,R}),
        // only coordinate r shifted. Classes outside nnz stay at zero, matching
        // the stationaryQlen of solver_amvald_forward.m rather than Qin, whose
        // rows for an empty class are not part of that point.
        // See _kb/06-solver-catalog.md (joint-dependence section).
        Matrix<T> jdterm(M, K, one);
        if (!jdscaling.empty()) {
            for (std::size_t r : nnz) {
                const bool closed = std::isfinite(Nc[r]);
                Matrix<T> arg(M, K, zero);
                for (std::size_t i = 0; i < M; ++i) {
                    for (std::size_t c : nnz) arg(i, c) = Qin(i, c);
                    T self = T(one + (closed ? selfArvl(i, r) : Qin(i, r)));
                    if (closed && linmethod)
                        self = T(self + num_traits<T>::from_double(Nc[r] - 1.0) * gamma[r](i, r));
                    arg(i, r) = self;
                }
                const std::vector<T> v = pfqn::pfqn_jdfun(arg, jdscaling, r);
                for (std::size_t i = 0; i < M; ++i) jdterm(i, r) = v[i];
            }
        }

        STeffOut = Matrix<T>(M, K, zero);
        for (std::size_t c : nnz)
            for (std::size_t i = 0; i < M; ++i)
                STeffOut(i, c) =
                    T(d.STchain(i, c) * lldterm(i, c) * msterm[i] * cdterm(i, c) * jdterm(i, c));

        // Wang-Sevcik queue line / fraction line: both REPLACE the arrival queue
        // length with an interpolation that needs STeff, so they run after it.
        if (method == "qli" || method == "fli") {
            const bool qli = (method == "qli");
            for (std::size_t i = 0; i < M; ++i) {
                const bool hol = L.stations[i].sched == SchedStrategy::HOL;
                for (std::size_t r : nnz) {
                    T tot = zero;
                    if (hol) {
                        for (std::size_t s2 : ehprio[r]) tot += Qin(i, s2);
                    } else {
                        for (std::size_t s2 : nnz) tot += Qin(i, s2);
                    }
                    if (Nc[r] == 1.0) {
                        totArvl(i, r) = T(tot - Qin(i, r));
                        continue;
                    }
                    const T num = T(STeffOut(i, r) * (one + tot - Qin(i, r)));
                    T den = zero;
                    for (std::size_t m = 0; m < M; ++m)
                        if (L.stations[m].sched == SchedStrategy::INF) den += STeffOut(m, r);
                    for (std::size_t m = 0; m < M; ++m) {
                        T totm = zero;
                        if (L.stations[m].sched == SchedStrategy::HOL) {
                            for (std::size_t s2 : ehprio[r]) totm += Qin(m, s2);
                        } else {
                            for (std::size_t s2 : nnz) totm += Qin(m, s2);
                        }
                        den += T(STeffOut(m, r) * (one + totm - Qin(m, r)));
                    }
                    if (den == zero) continue;
                    if (qli) {
                        const T f = T(one / num_traits<T>::from_double(Nc[r] - 1.0));
                        totArvl(i, r) = T(tot - f * (Qin(i, r) - num / den));
                    } else {
                        const T f = T(num_traits<T>::from_int(2) / num_traits<T>::from_double(Nc[r]));
                        totArvl(i, r) = T(tot - f * Qin(i, r) + num / den);
                    }
                }
            }
        }

        // Interlocked flow (Franks 1999, Eq. 4.7)
        // A request cannot queue behind work that its own submission caused, so the
        // arrival-instant queue drops the interlocked share of every other chain. The
        // own-class term is never removed. The matrix is empty for every model but the
        // layers of SolverLN, where it comes from the interlock path tables.
        if (!ILchain.empty()) {
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c : nnz) {
                    T ilq = zero;
                    for (std::size_t c2 : nnz)
                        if (c2 != c) ilq += ILchain(c, c2) * Qin(i, c2);
                    if (ilq > zero) {
                        T adj = T(totArvl(i, c) - ilq);
                        if (adj < selfArvl(i, c)) adj = selfArvl(i, c);
                        totArvl(i, c) = adj;
                    }
                }
        }

        Wout = Matrix<T>(M, K, zero);
        for (std::size_t c : nnz) {
            for (std::size_t i = 0; i < M; ++i) {
                const SchedStrategy sc = L.stations[i].sched;
                // Source residence-time rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                if (sc == SchedStrategy::EXT) continue;
                if (sc == SchedStrategy::INF) {
                    Wout(i, c) = STeffOut(i, c);
                    continue;
                }
                const double cs = L.stations[i].nservers;

                if (sc == SchedStrategy::PS) {
                    T corr = zero;
                    if (linmethod) {
                        for (std::size_t s : ccl)
                            corr += num_traits<T>::from_double(Nc[s]) * gamma[c](i, s);
                        corr -= gamma[c](i, c);
                    } else {
                        corr = T(num_traits<T>::from_double(Ntl - 1.0) * gamma[c](i, c));
                    }
                    if (opt.multiserver == "suri") {
                        // W = S (1 + Lm * suriFactor), against the arrival-instant
                        // queue seen by class c -- see the note below.
                        const T Lm = isopen[c] ? totArvlOpen(i, c)
                                               : T(totArvl(i, c) + corr);
                        Wout(i, c) = T(STeffOut(i, c) * (one + Lm * suri[i]));
                        continue;
                    }
                    if (opt.multiserver == "seidmann")
                        Wout(i, c) = T(STeffOut(i, c) * num_traits<T>::from_double(cs - 1.0));
                    if (isopen[c]) {
                        // PASTA no-discount rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                        Wout(i, c) = T(Wout(i, c) + STeffOut(i, c) * (one + totArvlOpen(i, c)));
                    } else {
                        // The arrival-instant queue is the one seen BY CLASS c:
                        // interpTot is a per-station total scaled by the GLOBAL
                        // (Ntot-1)/Ntot, so it cannot drop the arriving job's own
                        // chain and it charges a chain for jobs that never visit
                        // this station. totArvl(i,c) does drop it. Identical when
                        // the model has one chain -- see _kb/06-solver-catalog.md.
                        Wout(i, c) = T(Wout(i, c) + STeffOut(i, c) * (one + totArvl(i, c) + corr));
                    }
                    continue;
                }

                if (sc == SchedStrategy::DPS) {
                    // Discriminatory sharing: the arriving class waits behind
                    // each other class in proportion to the weight ratio.
                    const std::vector<T>& w = L.stations[i].schedparam;
                    if (w.size() != K)
                        throw UnsupportedError(
                            "solver_amvald: the DPS station '" + L.stations[i].name +
                            "' has no per-class weights");
                    T acc = zero;
                    if (cs > 1.0 && std::isfinite(cs))
                        acc = T(STeffOut(i, c) * num_traits<T>::from_double(cs - 1.0));
                    acc += T(STeffOut(i, c) * (one + selfArvl(i, c)));
                    for (std::size_t s : nnz) {
                        if (s == c) continue;
                        if (w[s] == w[c])
                            acc += T(STeffOut(i, c) * Qin(i, s));
                        else
                            acc += T(STeffOut(i, c) * Qin(i, s) * w[s] / w[c]);
                    }
                    Wout(i, c) = acc;
                    continue;
                }

                const bool fcfs_family = (sc == SchedStrategy::FCFS || sc == SchedStrategy::SIRO ||
                                          sc == SchedStrategy::LCFSPR);
                if (!fcfs_family && sc != SchedStrategy::HOL && sc != SchedStrategy::FCFSPRPRIO)
                    throw UnsupportedError(std::string("solver_amvald: ") + lang::sched_to_text(sc) +
                                           " scheduling is not implemented in this port");
                if (STeffOut(i, c) <= zero) continue;

                const std::vector<T>& Xarv = Xref.empty() ? Xin : Xref;

                // Preemptive-resume priority (PRIOMVA), Chandy-Lakshmi [ChaL83] applied
                // where it was derived: a job in service IS preempted by a higher-priority
                // arrival. Two terms separate this arm from the HOL path below.
                //  (a) no non-preemptive residual -- the lower-priority job found in
                //      service is preempted, so it delays nobody;
                //  (b) the tagged job's OWN service is interrupted too, so it is scaled
                //      by 1/(1-sigma_{k-1}) as well as the queued work.
                // Together: E[T_k] = E[S_k]/(1-sigma_{k-1}) + <queued>/(1-sigma_{k-1}).
                // Taken BEFORE the Bk / multiserver machinery because the reference
                // restricts the arm to single-server stations and refuses the rest.
                // Port of solver_amvald_forward.m (FCFSPRPRIO arm).
                if (sc == SchedStrategy::FCFSPRPRIO) {
                    if (cs > 1.0 && std::isfinite(cs))
                        throw UnsupportedError(
                            "solver_amvald: FCFSPRPRIO with more than one server. The "
                            "preemptive-resume priority arm (priomva) is implemented for "
                            "single-server stations only; use SolverCTMC or SolverSSA for "
                            "multiserver PRS.");

                    // higher-priority utilization seen at the arrival instant
                    T uh_prs = zero;
                    for (std::size_t h : hprio[c])
                        uh_prs += T(d.Vchain(i, h) * STeffOut(i, h) * T(Xarv[h] + tau[c][h]));
                    double vps = 1.0 - num_traits<T>::to_double(uh_prs);
                    vps = std::max(opt.tol, std::min(vps, 1.0 - opt.tol));
                    const T ps_prs = num_traits<T>::from_double(vps);

                    // work of EQUAL OR HIGHER priority already queued ahead
                    T queued = T(STeffOut(i, c) * (isopen[c] ? Qin(i, c) : selfArvl(i, c)));
                    for (std::size_t s : ehprio[c])
                        if (s != c) queued += T(STeffOut(i, s) * Qin(i, s));

                    Wout(i, c) = T((STeffOut(i, c) + queued) / ps_prs);
                    continue;
                }
                // Uchain_r rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                auto Ur = [&](std::size_t k, std::size_t s) -> T {
                    if (!(Xin[s] > zero)) return Uin(k, s);
                    return T(Uin(k, s) / Xin[s] * (Xarv[s] + tau[c][s]));
                };
                // prioScaling: the fraction of the server left by the strictly
                // higher-priority classes, clamped into [tol, 1-tol].
                auto prio_scaling = [&](std::size_t r) -> T {
                    if (sc != SchedStrategy::HOL) return one;
                    T uh = zero;
                    for (std::size_t h : hprio[r]) {
                        // shadow: Sevcik's server at the current population; cl: Eager-Lipscomb,
                        // the utilization seen at arrival, i.e. at population N - 1_r
                        const T xh = opt.np_priority == "shadow" ? Xarv[h] : T(Xarv[h] + tau[r][h]);
                        uh += T(d.Vchain(i, h) * STeffOut(i, h) * xh);
                    }
                    double v = 1.0 - num_traits<T>::to_double(uh);
                    v = std::max(opt.tol, std::min(v, 1.0 - opt.tol));
                    return num_traits<T>::from_double(v);
                };
                const T ps_r = prio_scaling(c);

                // Work of EQUAL OR HIGHER priority queued ahead of the arriving job, and the
                // residual of a strictly lower-priority job found in service, which HOL does not
                // preempt. This backlog is disjoint from the 1/ps_r inflation, which counts only
                // the higher-priority work that overtakes the job while it waits;
                // see _kb/06-solver-catalog.md
                auto hol_ehprio_backlog = [&](const std::vector<T>& Bkw) -> T {
                    T acc = zero;
                    for (std::size_t s : ehprio[c])
                        if (s != c) acc += T(STeffOut(i, s) * Qin(i, s) * Bkw[s]);
                    return acc;
                };
                auto hol_np_residual = [&](const std::vector<T>& Bkw) -> T {
                    if (opt.highvar == "hvmva") return zero; // hvmva already spans every class
                    T acc = zero;
                    for (std::size_t s : lprio[c])
                        acc += T(d.Vchain(i, s) * STeffOut(i, s) * (Xarv[s] + tau[c][s]) *
                                 STeffOut(i, s) * Bkw[s]);
                    return acc;
                };

                // Bk backlog-reduction rationale (FCFS/SIRO/LCFSPR vs HOL): see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                std::vector<T> Bk(K, one);
                if (cs > 1.0 && std::isfinite(cs)) {
                    T load = zero;
                    for (std::size_t s : nnz) {
                        const T dr = (sc == SchedStrategy::HOL) ? dcl[s] : ((s == c) ? dcl[c] : one);
                        load += dr * Xin[s] * d.Vchain(i, s) * STeffOut(i, s);
                    }
                    const bool light = num_traits<T>::to_double(load) < 0.75;
                    for (std::size_t s : nnz) {
                        const T dr = (sc == SchedStrategy::HOL) ? dcl[s] : ((s == c) ? dcl[c] : one);
                        T base = T(dr * Xin[s] * d.Vchain(i, s) * STeffOut(i, s));
                        if (sc == SchedStrategy::HOL && opt.multiserver != "softmin")
                            base = T(base / num_traits<T>::from_double(cs));
                        if (light) {
                            Bk[s] = base;
                        } else {
                            const unsigned e = static_cast<unsigned>(
                                std::llround(cs) - (sc == SchedStrategy::HOL ? 0 : 1));
                            Bk[s] = num_pow_int(base, e);
                        }
                    }
                }

                // single-server-with-scaling branch rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                if (cs == 1.0) {
                    T w = zero;
                    if (opt.highvar == "hvmva") {
                        T usum = zero;
                        for (std::size_t s : ccl) usum += Ur(i, s);
                        w = T(STeffOut(i, c) * (one - usum));
                        for (std::size_t s : ccl)
                            w += T(STeffOut(i, s) / prio_scaling(s) * Ur(i, s) *
                                   (one + d.SCVchain(i, s)) / num_traits<T>::from_int(2));
                    } else if (opt.highvar != "default") {
                        throw UnsupportedError("solver_amvald: highvar '" + opt.highvar +
                                               "' is not implemented");
                    } else {
                        w = STeffOut(i, c);
                    }
                    if (sc == SchedStrategy::HOL) {
                        // priority-branch charging rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                        const std::vector<T> ones(K, one);
                        w += T((STeffOut(i, c) * (isopen[c] ? Qin(i, c) : selfArvl(i, c)) +
                                hol_ehprio_backlog(ones) + hol_np_residual(ones)) / ps_r);
                    } else {
                        w += STeffOut(i, c) * (isopen[c] ? Qin(i, c) : selfArvl(i, c));
                        for (std::size_t s : nnz)
                            if (s != c) w += STeffOut(i, s) * Qin(i, s);
                    }
                    Wout(i, c) = w;
                    continue;
                }

                // Multiserver.
                if (opt.multiserver == "suri") {
                    T Lm = isopen[c] ? T(deltaclass[c] * Qin(i, c)) : selfArvl(i, c);
                    if (sc != SchedStrategy::HOL)
                        for (std::size_t s : nnz)
                            if (s != c) Lm += Qin(i, s);
                    if (sc == SchedStrategy::HOL) {
                        Lm = isopen[c] ? Qin(i, c) : selfArvl(i, c);
                        const std::vector<T> ones(K, one);
                        Wout(i, c) = T(STeffOut(i, c) +
                                       ((STeffOut(i, c) * Lm + hol_ehprio_backlog(ones)) * suri[i] +
                                        hol_np_residual(ones)) / ps_r);
                        continue;
                    }
                    Wout(i, c) = T(STeffOut(i, c) / ps_r + STeffOut(i, c) * Lm * suri[i] / ps_r);
                    continue;
                }
                if (opt.multiserver == "softmin") {
                    T w = STeffOut(i, c);
                    if (sc == SchedStrategy::HOL)
                        w += T((STeffOut(i, c) * (isopen[c] ? Qin(i, c) : selfArvl(i, c)) * Bk[c] +
                                hol_ehprio_backlog(Bk) + hol_np_residual(Bk)) / ps_r);
                    else
                        w += T(STeffOut(i, c) * (isopen[c] ? Qin(i, c) : selfArvl(i, c)) * Bk[c] / ps_r);
                    if (sc != SchedStrategy::HOL)
                        for (std::size_t s : nnz)
                            if (s != c) w += STeffOut(i, s) * Bk[s] * Qin(i, s);
                    Wout(i, c) = w;
                    continue;
                }
                // default / seidmann
                T w = T(STeffOut(i, c) * num_traits<T>::from_double(cs - 1.0));
                w += STeffOut(i, c);
                if (sc == SchedStrategy::HOL) {
                    w += T((STeffOut(i, c) * (isopen[c] ? Qin(i, c) : selfArvl(i, c)) * Bk[c] +
                            hol_ehprio_backlog(Bk) + hol_np_residual(Bk)) / ps_r);
                } else {
                    w += STeffOut(i, c) *
                         (isopen[c] ? T(deltaclass[c] * Qin(i, c)) : selfArvl(i, c)) * Bk[c];
                    for (std::size_t s : nnz)
                        if (s != c) w += STeffOut(i, s) * Bk[s] * Qin(i, s);
                }
                Wout(i, c) = w;
            }
        }
    };

    // one inner fixed point at a given population vector
    auto inner_loop = [&](Matrix<T>& Q, std::vector<T>& X, Matrix<T>& U,
                          const std::vector<double>& Nc, Matrix<T>& Wout, Matrix<T>& STeffOut) {
        int iter = 0;
        Matrix<T> Qprev = Q;
        while (true) {
            Qprev = Q;
            const std::vector<T> Xprev = X;
            const Matrix<T> Uprev = U;
            ++iter;
            forward(Qprev, Xprev, Uprev, Nc, Wout, STeffOut);
            ++totiter;
            if (line::util::LineConsole::owns_log()) {
                double resid = 0.0;
                for (std::size_t i = 0; i < Q.rows(); ++i)
                    for (std::size_t c = 0; c < Q.cols(); ++c)
                        resid = std::max(resid, std::abs(num_traits<T>::to_double(Q(i, c)) -
                                                         num_traits<T>::to_double(Qprev(i, c))));
                line::util::LineConsole::iter(totiter,
                                              "AMVA sweep %d: queue-length residual %.3e",
                                              totiter, resid);
            }
            if (totiter >= max_totiter) break;
            for (std::size_t c : nnz) {
                T sw = zero;
                for (std::size_t i = 0; i < M; ++i) sw += Wout(i, c);
                if (sw == zero) {
                    X[c] = zero;
                } else if (!isopen[c]) {
                    T cyc = zero;
                    for (std::size_t i = 0; i < M; ++i) cyc += d.Vchain(i, c) * Wout(i, c);
                    if (cyc != zero)
                        X[c] = T(omicron * num_traits<T>::from_double(Nc[c]) / cyc +
                                 (one - omicron) * Xprev[c]);
                }
                // an open chain's X is its arrival rate and stays put; only its
                // queue lengths and utilizations follow from the residence times
                for (std::size_t i = 0; i < M; ++i) {
                    Q(i, c) = T(omicron * X[c] * d.Vchain(i, c) * Wout(i, c) +
                                (one - omicron) * Qprev(i, c));
                    U(i, c) = T(omicron * d.Vchain(i, c) * STeffOut(i, c) * X[c] +
                                (one - omicron) * Uprev(i, c));
                }
            }
            if (iter >= 2) {
                double err = 0.0;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t c = 0; c < K; ++c)
                        err = std::max(err, std::fabs(num_traits<T>::to_double(Q(i, c) - Qprev(i, c))));
                if (err <= opt.iter_tol) break;
            }
            if (iter > inner_cap) break;
        }
        return Qprev;
    };

    Matrix<T> QouterPrev = Qchain;
    int outer = 0;
    converged = false;
    while (true) {
        ++outer;
        QouterPrev = Qchain;
        const std::vector<T> XouterPrev = Xchain;
        // baseline the tau differences below are taken against; empty when no recursion runs
        if (linmethod) Xref = XouterPrev;

        // Linearizer correction: solve at each reduced population N - e_s
        if (linmethod && std::isfinite(Nt) && Nt > 0.0) {
            for (std::size_t s = 0; s < K; ++s) {
                if (!std::isfinite(d.Nchain[s])) continue;
                std::vector<double> Ns = d.Nchain;
                Ns[s] -= 1.0;
                const T shrink = num_traits<T>::from_double((Nt - 1.0) / Nt);
                Matrix<T> Qs = Qchain;
                std::vector<T> Xs = Xchain;
                Matrix<T> Us = Uchain;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t c = 0; c < K; ++c) {
                        Qs(i, c) = T(Qs(i, c) * shrink);
                        Us(i, c) = T(Us(i, c) * shrink);
                    }
                for (std::size_t c = 0; c < K; ++c) Xs[c] = T(Xs[c] * shrink);
                Matrix<T> Ws, STs;
                const std::vector<T> Xs_in = Xs;
                const Matrix<T> Qs_prev = inner_loop(Qs, Xs, Us, Ns, Ws, STs);
                // tau(s,r) rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
                for (std::size_t c : nnz) tau[s][c] = T(Xs[c] - XouterPrev[c]);
                (void)Xs_in;
                // ONLY 'lin' TAKES THE PER-CLASS CORRECTION. 'qdlin' takes the
                // class-AGGREGATE one and stores it in column 0, leaving the
                // rest of the row zero, which is what the reference does:
                // solver_amvald.m allocates the (K,M,K) per-class array for
                // qdlin but writes `gamma(s,k) = sum_r Q_s(k,r)/(Nt-1) -
                // sum_r Q(k,r)/Nt` into it with two subscripts, and MATLAB
                // linear-indexes that to (s,k,1). Every reader below still
                // indexes gamma PER CLASS, so the correction that reaches the
                // residence time is N_0*gamma(r,i,0) - [r==0]*gamma(r,i,0). It
                // coincides with the queue-dependent AMVA form (Nt-1)*gamma_agg
                // iff K = 1, so a single-chain model is unaffected and lin and
                // qdlin stay bit-identical there. This port filled it per class
                // for both methods until 2026-09-04, which made C++ qdlin an
                // alias of C++ lin and put it at odds with MATLAB, the JAR and
                // native Python on every multichain model. Do not "restore" it
                // without re-baselining the AMVA goldens in all four codebases;
                // see _kb/06-solver-catalog.md.
                if (method == "qdlin") {
                    for (std::size_t i = 0; i < M; ++i) {
                        T qs = zero, qo = zero;
                        for (std::size_t c = 0; c < K; ++c) {
                            qs += Qs_prev(i, c);
                            qo += QouterPrev(i, c);
                        }
                        // Nt = 1 leaves the reduced population empty and the
                        // reference divides by zero there; native Python guards
                        // it with a zero and this port follows Python.
                        gamma[s](i, 0) =
                            (Nt > 1.0) ? T(qs / num_traits<T>::from_double(Nt - 1.0) -
                                           qo / num_traits<T>::from_double(Nt))
                                       : zero;
                    }
                } else {
                    for (std::size_t i = 0; i < M; ++i)
                        for (std::size_t c : nnz)
                            if (std::isfinite(d.Nchain[c]) && Ns[c] > 0.0)
                                gamma[s](i, c) =
                                    T(Qs_prev(i, c) / num_traits<T>::from_double(Ns[c]) -
                                      QouterPrev(i, c) / num_traits<T>::from_double(d.Nchain[c]));
                }
                if (totiter >= max_totiter) break;
            }
        }
        if (totiter >= max_totiter) break;

        Matrix<T> Wtmp, STtmp;
        const Matrix<T> Qprev = inner_loop(Qchain, Xchain, Uchain, d.Nchain, Wtmp, STtmp);
        Wchain = Wtmp;
        STeff = STtmp;

        double err = 0.0;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < K; ++c)
                err = std::max(err, std::fabs(num_traits<T>::to_double(Qchain(i, c) - QouterPrev(i, c))));
        converged = err <= opt.iter_tol;
        if (outer >= 2 && converged) break;
        if (outer >= inner_cap || totiter > max_totiter) break;
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));

    // renormalise utilizations past one, as the reference does
    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy sc = L.stations[i].sched;
        if (!(sc == SchedStrategy::FCFS || sc == SchedStrategy::SIRO || sc == SchedStrategy::PS ||
              sc == SchedStrategy::LCFSPR || sc == SchedStrategy::DPS || sc == SchedStrategy::HOL))
            continue;
        T usum = zero;
        for (std::size_t c = 0; c < K; ++c) usum += Uchain(i, c);
        if (num_traits<T>::to_double(usum) <= 1.0) continue;
        T den = zero;
        for (std::size_t c = 0; c < K; ++c) den += d.Vchain(i, c) * STeff(i, c) * Xchain[c];
        if (den == zero) continue;
        for (std::size_t c = 0; c < K; ++c) {
            if (!(d.Vchain(i, c) * STeff(i, c) > zero)) continue;
            Uchain(i, c) = T(one * d.Vchain(i, c) * STeff(i, c) * Xchain[c] / den);
        }
    }

    Matrix<T> Rchain(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c)
            if (Tchain(i, c) != zero) Rchain(i, c) = T(Qchain(i, c) / Tchain(i, c));
    for (std::size_t c = 0; c < K; ++c) {
        if (d.Nchain[c] != 0.0) continue;
        Xchain[c] = zero;
        for (std::size_t i = 0; i < M; ++i) {
            Uchain(i, c) = zero;
            Rchain(i, c) = zero;
            Tchain(i, c) = zero;
        }
    }

    const ClassResults<T> cr =
        sn_deaggregate_chain_results(L, d, Matrix<T>(), Matrix<T>(), Rchain, Tchain, Xchain);
    MvaSolution<T> out;
    out.Q = cr.Q;
    out.U = cr.U;
    out.R = cr.R;
    out.Tp = cr.Tp;
    out.C = cr.C;
    out.X = cr.X;

    // A station with limited class or joint dependence reports utilization as
    // T*S/peak from the declared peak rate scaling, matching the T*S/c
    // convention of an ordinary multiserver station (solver_amvald.m:255-291).
    //
    // THE PEAK IS A MODEL INPUT, NOT A DEFAULT. Without it this pass has no
    // normalizer, and returning the unnormalized column -- or, as it did before,
    // a column of zeros -- is worse than refusing: both read as a utilization.
    // SolverCTMC, solver_nc_conv and both SSA engines already refuse it by name
    // here, and MATLAB's getLimitedClassDependencePeak.m refuses it in the
    // model. This check sits ahead of the closed-model guard below because a
    // missing peak is a defect whether or not the pass would have run.
    for (std::size_t ist = 0; ist < M; ++ist) {
        if (L.stations[ist].cdscaling && L.stations[ist].cdscalingpeak.empty())
            throw InputError(
                "SolverMVA: station '" + L.stations[ist].name +
                "' declares class-dependent service without a peak rate. Utilization at a "
                "class-dependent station is reported as T*S/peak, so pass the peak to "
                "setClassDependence");
        if (L.stations[ist].jdscaling && L.stations[ist].jdscalingpeak.empty())
            throw InputError(
                "SolverMVA: station '" + L.stations[ist].name +
                "' declares joint-dependent service without a peak rate; pass the peak to "
                "setJointDependence");
    }
    bool anyOpen = false;
    for (const auto& c : L.classes)
        if (std::isinf(c.population)) anyOpen = true;
    if (!anyOpen) {
        const std::size_t Kcls = L.nclasses;
        for (std::size_t ist = 0; ist < M; ++ist) {
            const std::vector<T>* peak = nullptr;
            if (L.stations[ist].cdscaling)
                peak = &L.stations[ist].cdscalingpeak;
            else if (L.stations[ist].jdscaling)
                peak = &L.stations[ist].jdscalingpeak;
            if (peak == nullptr) continue;
            for (std::size_t k = 0; k < Kcls; ++k) {
                const double rate = num_traits<T>::to_double(L.rates(ist, k));
                const T bmax = k < peak->size() ? (*peak)[k] : zero;
                if (std::isfinite(rate) && rate > 0.0 && bmax > zero)
                    out.U(ist, k) = T(out.Tp(ist, k) / L.rates(ist, k) / bmax);
                else
                    out.U(ist, k) = zero;
            }
        }
    }
    out.method = method;
    out.iter = totiter;
    // through DETAIL, not STEP: SolverLN runs this analyzer once per layer per
    // iteration, inside its own run, so a step line here repeats hundreds of
    // times; detail() collapses repeats of the same shape
    {
        char buf[96];
        std::snprintf(buf, sizeof(buf), "AMVA finished after %d sweeps", totiter);
        line::util::LineConsole::detail(buf);
    }
    return out;
}

// ---------------------------------------------------------------------------
// solver_amva: method resolution and the product-form fast paths
// ---------------------------------------------------------------------------

/**
 * The @c amva.* spellings of solver_amva.m, mapped onto the bare method names.
 *
 * The reference accepts both, and a model or a test may use either, so the
 * alias table is part of the interface rather than a convenience: an
 * unrecognised name must reach the method switch and be refused there by its
 * own name, not silently resolved to `default`.
 */
inline std::string amva_method_alias(const std::string& m) {
    if (m == "amva.qli") return "qli";
    if (m == "amva.qd" || m == "amva.qdamva" || m == "qdamva") return "qd";
    if (m == "amva.aql") return "aql";
    if (m == "amva.qsa") return "qsa";
    if (m == "amva.qdaql") return "qdaql";
    if (m == "amva.tay") return "tay";
    if (m == "amva.scat") return "scat";
    if (m == "amva.lin") return "lin";
    if (m == "amva.qdlin") return "qdlin";
    if (m == "amva.fli") return "fli";
    if (m == "amva.bs") return "bs";
    if (m == "amva.ab") return "ab";
    if (m == "amva.schmidt") return "schmidt";
    if (m == "amva.schmidt-ext") return "schmidt-ext";
    if (m == "amva.lcp") return "lcp";
    if (m == "amva.chow") return "chow";
    if (m == "amva.pamb") return "pamb";
    if (m == "amva.pami") return "pami";
    if (m == "amva.pamt") return "pamt";
    if (m == "amva.clust") return "clust";
    if (m == "amva.dmlin") return "dmlin";
    if (m == "amva.priomva") return "priomva";
    if (m == "amva.marie") return "marie";
    return m;
}

template <class T>
MvaSolution<T> solver_amva(const qn::NetworkStruct<T>& L, const ChainDemands<T>& d, MvaOptions opt,
                           const Matrix<T>& init_sol, bool& converged) {
    using qn::SchedStrategy;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, C = L.nchains;
    opt.iter_max = std::min(opt.iter_max, 10000);
    converged = true;

    std::string method = amva_method_alias(opt.method);
    if (method == "default" || method == "amva") {
        double Nsum = 0.0;
        bool anysmall = false;
        for (std::size_t c = 0; c < C; ++c) {
            Nsum += d.Nchain[c];
            if (d.Nchain[c] < 1.0) anysmall = true;
        }
        if (Nsum <= 2.0 || anysmall) {
            method = "qd";
        } else {
            bool anyfinite = false, allone = true;
            for (const auto& s : L.stations)
                if (std::isfinite(s.nservers)) {
                    anyfinite = true;
                    if (s.nservers != 1.0) allone = false;
                }
            // empty finite-server max rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
            method = (anyfinite && allone) ? "egflin" : "lin";
        }
    }

    // The closed-population AMVA family (Bard-Schweitzer, SQNI, Tay, SCAT, AQL, QSA,
    // Bard LCP, Chow SA, Hsieh-Lam PAM, clustering, Improved Linearizer,
    // Akyildiz-Bolch, Schmidt) lives ONLY in the product-form branch below. The same
    // predicate the report gates on decides here, so a name the report offers is a
    // name that runs and a name it withholds errors rather than falling through to
    // solver_amvald and returning the qd-family answer under a method the caller did
    // not ask for.
    {
        const std::string amva_reason = mva_closed_population_reason(L, method);
        if (!amva_reason.empty()) throw UnsupportedError(amva_reason);
    }

    // MATLAB's "trivial models" early exit tests sn_has_homogeneous_scheduling,
    // which reduces to nstations == 1; see the note on that predicate.
    if (L.has_homogeneous_scheduling(SchedStrategy::INF)) {
        MvaOptions o2 = opt;
        o2.multiserver = "default";
        return solver_amvald(L, d, o2, method, converged, init_sol);
    }

    // ab / schmidt / schmidt-ext ARE the class-dependent FCFS algorithms and live only in
    // the product-form branch below, so the het-FCFS exclusion must not divert them.
    const bool het_fcfs_own =
        (method == "ab" || method == "schmidt" || method == "schmidt-ext");
    const bool cond1 = L.has_product_form_not_het_fcfs() ||
                       (het_fcfs_own && L.has_product_form_not_het_fcfs(false));
    bool cond2 = true;  // ~sn_has_load_dependence: the pf kernels ignore lldscaling
    for (const auto& st : L.stations)
        if (!st.lldscaling.empty()) cond2 = false;
    // MATLAB keeps ONE mixed model in this branch (solver_amva.m:146): a strict
    // product-form network under 'lin', where pfqn_linearizermx solves the open classes
    // in closed form and inflates the closed demands by their utilization. Every other
    // open model goes to solver_amvald, the only path carrying the open corrections.
    const bool cond3 = !L.has_open_classes() || (L.has_product_form() && method == "lin");
    if (!(cond1 && cond2 && cond3)) return solver_amvald(L, d, opt, method, converged, init_sol);

    // Interlocked flow (Franks 1999, Eq. 4.7): only solver_amvald carries the correction, so
    // an interlocked model goes there rather than to the product-form kernels below, which
    // have no interlock term and would drop it silently. Mirrors solver_amva.m and the JAR.
    if (!sn_interlock_chain<T>(L, opt.interlock).empty())
        return solver_amvald(L, d, opt, method, converged, init_sol);

    const PfChainParams<T> pf = sn_get_product_form_chain_params(L, d);
    if (pf.queue_stations.empty()) return solver_amvald(L, d, opt, method, converged, init_sol);

    const std::size_t nq = pf.queue_stations.size(), nz = pf.delay_stations.size();
    // An open chain has no population: it takes pfqn's kOpenClass sentinel, which stands
    // in for MATLAB's Inf, and llround(Inf) is undefined behaviour rather than a marker.
    std::vector<int> N(C, 0);
    for (std::size_t c = 0; c < C; ++c)
        N[c] = std::isfinite(d.Nchain[c]) ? static_cast<int>(std::llround(d.Nchain[c]))
                                          : pfqn::kOpenClass;
    // Nt is read only by the closed-only kernels, which no open model reaches, and an
    // exact field has no value for from_double(Inf); an open chain is left at zero.
    std::vector<T> Nt(C, zero);
    for (std::size_t c = 0; c < C; ++c)
        if (std::isfinite(d.Nchain[c])) Nt[c] = num_traits<T>::from_double(d.Nchain[c]);

    // The scheduling of the queueing stations, in the api layer's enum.
    std::vector<pfqn::SchedStrategy> types;
    for (std::size_t a = 0; a < nq; ++a) {
        switch (L.stations[pf.queue_stations[a] - 1].sched) {
            case SchedStrategy::PS: types.push_back(pfqn::SchedStrategy::PS); break;
            case SchedStrategy::INF: types.push_back(pfqn::SchedStrategy::INF); break;
            default: types.push_back(pfqn::SchedStrategy::FCFS); break;
        }
    }

    // Warm start, aligned to the queueing-station rows; a malformed seed is
    // dropped rather than repaired, as the reference does.
    Matrix<T> Q0;
    if (init_sol.rows() == L.nstations && init_sol.cols() == C) {
        Q0 = Matrix<T>(nq, C, zero);
        bool bad = false;
        for (std::size_t a = 0; a < nq && !bad; ++a)
            for (std::size_t c = 0; c < C; ++c) {
                Q0(a, c) = init_sol(pf.queue_stations[a] - 1, c);
                if (Q0(a, c) < zero) {
                    bad = true;
                    break;
                }
            }
        if (bad) Q0 = Matrix<T>();
    }

    // Seidmann-scaling exemption rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    const bool direct_ms = (method == "ab" || method == "schmidt" || method == "schmidt-ext");
    const bool seidmann = (opt.multiserver == "default" || opt.multiserver == "seidmann");

    Matrix<T> Dm = pf.D, Zm = pf.Z;
    if (!direct_ms) {
        if (seidmann) {
            // move the (c-1)/c share of each multiserver demand into the first
            // delay station
            for (std::size_t a = 0; a < nq; ++a) {
                const double c = pf.S[a];
                if (!std::isfinite(c)) continue;
                for (std::size_t k = 0; k < C; ++k)
                    Dm(a, k) = T(Dm(a, k) / num_traits<T>::from_double(c));
                if (Zm.rows() > 0)
                    for (std::size_t k = 0; k < C; ++k)
                        Zm(0, k) =
                            T(Zm(0, k) + pf.D(a, k) * num_traits<T>::from_double((c - 1.0) / c));
            }
        } else if (opt.multiserver == "softmin") {
            return solver_amvald(L, d, opt, method, converged, init_sol);
        }
        // conway, krzesinski and erlang leave the demands alone: they carry the
        // multiplicity into the algorithm itself.
    }

    // Think time per class, which is what the single-server api routines take.
    std::vector<T> Zsum(C, zero);
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t a = 0; a < Zm.rows(); ++a) Zsum[c] += Zm(a, c);

    bool allone = true, anyinf = false;
    for (double s : pf.S) {
        if (!std::isfinite(s)) anyinf = true;
        else if (s != 1.0) allone = false;
    }

    // The per-method solve fills these, over the QUEUEING stations only.
    Matrix<T> Qq(nq, C, zero), Uq(nq, C, zero), Qz(nz, C, zero);
    std::vector<T> X(C, zero);
    int iters = 0;
    bool have_delay_rows = false;  // set by the methods that solve the delays too

    if (method == "sqni") {
        // square-root approximation gate rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
        if (!(M == 2 && nq == 1 && nz == 1))
            throw UnsupportedError(
                "solver_amva: method 'sqni' applies only to a model of one queueing station and "
                "one infinite server");
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "solver_amva: method 'sqni' solves a quadratic and evaluates a square root, "
                "which exact arithmetic has no representation for; use the double or real "
                "backend");
        } else {
            std::vector<T> Lv(C, zero);
            for (std::size_t c = 0; c < C; ++c) Lv[c] = Dm(0, c);
            const pfqn::SqniResult<T> r = pfqn::pfqn_sqni(Nt, Lv, Zsum);
            for (std::size_t c = 0; c < C; ++c) {
                Qq(0, c) = r.Q[c];
                Uq(0, c) = r.U[c];
                X[c] = r.X[c];
            }
            iters = 1;
        }
    } else if (method == "bs") {
        std::vector<pfqn::AmvaSched> bstype;
        for (std::size_t a = 0; a < nq; ++a)
            bstype.push_back(L.stations[pf.queue_stations[a] - 1].sched == SchedStrategy::PS
                                 ? pfqn::AmvaSched::PS
                                 : pfqn::AmvaSched::FCFS);
        const pfqn::AmvaResult<T> r =
            pfqn::pfqn_bs(Dm, Nt, Zsum, bstype, opt.tol, static_cast<std::size_t>(opt.iter_max), Q0);
        Qq = r.QN;
        Uq = r.UN;
        X = r.XN;
        iters = static_cast<int>(r.iterations);
    } else if (method == "lcp" || method == "chow") {
        // Bard LCP and the Chow Second Approximation built on it: both are
        // Bard-Schweitzer variants in the arrival-instant estimate, so they take
        // the same Seidmann treatment and the same scheduling tags as 'bs'.
        std::vector<pfqn::AmvaSched> lctype;
        for (std::size_t a = 0; a < nq; ++a)
            lctype.push_back(L.stations[pf.queue_stations[a] - 1].sched == SchedStrategy::PS
                                 ? pfqn::AmvaSched::PS
                                 : pfqn::AmvaSched::FCFS);
        const pfqn::AmvaResult<T> r =
            (method == "lcp")
                ? pfqn::pfqn_lcp(Dm, Nt, Zsum, lctype, opt.tol,
                                 static_cast<std::size_t>(opt.iter_max), Q0)
                : pfqn::pfqn_chow(Dm, Nt, Zsum, lctype, opt.tol,
                                  static_cast<std::size_t>(opt.iter_max), Q0);
        Qq = r.QN;
        Uq = r.UN;
        X = r.XN;
        iters = static_cast<int>(r.iterations);
    } else if (method == "pamb" || method == "pami" || method == "pamt") {
        // Hsieh-Lam proportional approximations, noniterative
        const pfqn::PamVariant pv = (method == "pamb")   ? pfqn::PamVariant::Basic
                                    : (method == "pami") ? pfqn::PamVariant::Improved
                                                         : pfqn::PamVariant::Two;
        const pfqn::AmvaResult<T> r = pfqn::pfqn_pam(Dm, Nt, Zsum, pv);
        Qq = r.QN;
        Uq = r.UN;
        X = r.XN;
        iters = 1;
    } else if (method == "clust") {
        // de Souza e Silva-Lavenberg-Muntz clustering approximation, with the
        // decomposition derived automatically from the PAMB utilizations.
        const pfqn::AmvaResult<T> r = pfqn::pfqn_clust(
            Dm, Nt, Zsum, std::vector<std::vector<std::size_t> >(),
            std::vector<std::vector<std::size_t> >(), pfqn::ClustInner::Linearizer, opt.tol,
            static_cast<std::size_t>(opt.iter_max));
        Qq = r.QN;
        Uq = r.UN;
        X = r.XN;
        iters = static_cast<int>(r.iterations);
    } else if (method == "dmlin") {
        // de Souza e Silva-Muntz Improved Linearizer: the Linearizer fixed point
        // with the Delta-terms pre-aggregated, so it agrees with 'lin' by
        // construction and only the cost differs.
        const pfqn::LinearizerResult<T> r =
            pfqn::pfqn_dmlin(Dm, N, Zm, types, opt.tol, opt.iter_max, Q0);
        Qq = r.Q;
        Uq = r.U;
        X = r.X;
        iters = r.totiter;
    } else if (method == "tay") {
        // Between 'bs' and 'aql', which is the reference's own branch order in
        // solver_amva.m. MATLAB raises here rather than approximating: the
        // elasticity equations are derived for single servers throughout, and
        // there is no multiserver correction to fall back on.
        if (L.has_multi_server())
            throw UnsupportedError(
                "solver_amva: Tay's approximation is defined for single-server stations; use "
                "'default' or 'lin'");
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "solver_amva: method 'tay' stops on a tolerance evaluated with transcendentals, "
                "which exact arithmetic has no representation for; use the double or real "
                "backend");
        } else {
            const pfqn::AmvaResult<T> r = pfqn::pfqn_tay(
                Dm, Nt, Zsum, opt.tol, static_cast<std::size_t>(opt.iter_max), Q0);
            Qq = r.QN;
            Uq = r.UN;
            X = r.XN;
            iters = static_cast<int>(r.iterations);
        }
    } else if (method == "scat") {
        // Neuse-Chandy SCAT: the Linearizer fixed point with one Delta refresh
        // instead of three. Unlike tay/aql/qsa this takes multiserver stations,
        // because they reach it already Seidmann-scaled, exactly as for 'bs'.
        const pfqn::LinearizerResult<T> r =
            pfqn::pfqn_scat(Dm, N, Zm, types, opt.tol, opt.iter_max, Q0);
        Qq = r.Q;
        Uq = r.U;
        X = r.X;
        iters = r.totiter;
    } else if (method == "aql") {
        // MATLAB raises here rather than approximating: the AQL recursion has
        // no multiserver correction at all.
        if (L.has_multi_server())
            throw UnsupportedError(
                "solver_amva: AQL cannot handle multi-server stations; use 'default' or 'lin'");
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "solver_amva: method 'aql' stops on a relative tolerance evaluated with "
                "transcendentals, which exact arithmetic has no representation for; use the "
                "double or real backend");
        } else {
            const pfqn::AmvaResult<T> r =
                pfqn::pfqn_aql(Dm, Nt, Zsum, opt.tol, static_cast<std::size_t>(opt.iter_max));
            Qq = r.QN;
            Uq = r.UN;
            X = r.XN;
            iters = static_cast<int>(r.iterations);
        }
    } else if (method == "qsa") {
        // MATLAB raises here rather than approximating: the QSA core carries an
        // aggregate queue length with no multiserver correction at all.
        if (L.has_multi_server())
            throw UnsupportedError(
                "solver_amva: QSA cannot handle multi-server stations; use 'default' or 'lin'");
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "solver_amva: method 'qsa' stops on a residual tolerance evaluated with "
                "transcendentals, which exact arithmetic has no representation for; use the "
                "double or real backend");
        } else {
            // Dm holds the queueing stations only, the delays already folded
            // into Zsum, so every row here is a queueing centre.
            const pfqn::AmvaResult<T> r = pfqn::pfqn_qsa(
                Dm, Nt, Zsum, std::vector<pfqn::AmvaSched>(nq, pfqn::AmvaSched::PS), opt.tol,
                static_cast<std::size_t>(opt.iter_max));
            Qq = r.QN;
            Uq = r.UN;
            X = r.XN;
            iters = static_cast<int>(r.iterations);
        }
    } else if (direct_ms) {
        // The demands of the delays are stacked ON TOP of the queues, which is
        // the [Z0; L0] layout these three routines expect.
        Matrix<T> Dfull(nz + nq, C, zero), Vfull(nz + nq, C, num_traits<T>::from_int(1));
        std::vector<int> nsfull;
        std::vector<pfqn::SchedStrategy> schedfull;
        Matrix<int> Sfull(nz + nq, 1, 1);
        for (std::size_t a = 0; a < nz; ++a) {
            for (std::size_t c = 0; c < C; ++c) Dfull(a, c) = pf.Z(a, c);
            nsfull.push_back(1);  // an infinite server is marked by its discipline
            schedfull.push_back(pfqn::SchedStrategy::INF);
            Sfull(a, 0) = 1;
        }
        for (std::size_t a = 0; a < nq; ++a) {
            for (std::size_t c = 0; c < C; ++c) {
                Dfull(nz + a, c) = pf.D(a, c);
                Vfull(nz + a, c) = d.Vchain(pf.queue_stations[a] - 1, c) > zero
                                       ? num_traits<T>::from_int(1)
                                       : zero;
            }
            const int c_i = std::isfinite(pf.S[a]) ? static_cast<int>(std::llround(pf.S[a])) : 1;
            nsfull.push_back(c_i);
            Sfull(nz + a, 0) = c_i;
            schedfull.push_back(types[a]);
        }
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "solver_amva: the Akyildiz-Bolch and Schmidt multiserver methods use marginal "
                "weights with non-integer powers and floors, which exact arithmetic has no "
                "representation for; use the double or real backend");
        } else if (method == "ab") {
            const pfqn::AbAmvaResult<T> r = pfqn::pfqn_ab_amva(
                Dfull, N, Vfull, nsfull, schedfull, false, pfqn::AbMarginalMethod::Ab);
            for (std::size_t a = 0; a < nq; ++a)
                for (std::size_t c = 0; c < C; ++c) {
                    Qq(a, c) = r.QN(nz + a, c);
                    Uq(a, c) = r.UN(nz + a, c);
                }
            for (std::size_t a = 0; a < nz; ++a)
                for (std::size_t c = 0; c < C; ++c) Qz(a, c) = r.QN(a, c);
            X = r.XN;
            iters = static_cast<int>(r.totiter);
        } else if (method == "schmidt") {
            const pfqn::SchmidtResult<T> r = pfqn::pfqn_schmidt(Dfull, N, Sfull, schedfull, Vfull);
            for (std::size_t a = 0; a < nq; ++a)
                for (std::size_t c = 0; c < C; ++c) Qq(a, c) = r.QN(nz + a, c);
            for (std::size_t a = 0; a < nz; ++a)
                for (std::size_t c = 0; c < C; ++c) Qz(a, c) = r.QN(a, c);
            X = r.XN;
            iters = 1;
        } else {
            // One predicate for the gate and the run, asked about the numbers THIS
            // arm passes: pfqn_schmidt_ext forms its alpha correction from the
            // network with one class-r customer tagged, and a chain holding no
            // customer has none to tag.
            {
                std::vector<double> sx_n(C, 0.0);
                for (std::size_t c = 0; c < C; ++c) sx_n[c] = d.Nchain[c];
                std::vector<bool> sx_fcfs;
                for (std::size_t a = 0; a < nq; ++a)
                    sx_fcfs.push_back(L.stations[pf.queue_stations[a] - 1].sched ==
                                      SchedStrategy::FCFS);
                const std::string sx_reason =
                    mva_schmidt_ext_reason(sx_n, sx_fcfs, "schmidt-ext");
                if (!sx_reason.empty()) throw UnsupportedError(sx_reason);
            }
            const pfqn::SchmidtExtResult<T> r = pfqn::pfqn_schmidt_ext(Dfull, N, Sfull, schedfull);
            for (std::size_t a = 0; a < nq; ++a)
                for (std::size_t c = 0; c < C; ++c) Qq(a, c) = r.QN(nz + a, c);
            for (std::size_t a = 0; a < nz; ++a)
                for (std::size_t c = 0; c < C; ++c) Qz(a, c) = r.QN(a, c);
            X = r.XN;
            iters = 1;
        }
        have_delay_rows = true;
        // U = X D / c, which is what the reference recomputes for both Schmidt
        // variants rather than reading it back from the algorithm.
        for (std::size_t a = 0; a < nq; ++a) {
            const double c_i = std::isfinite(pf.S[a]) ? pf.S[a] : 1.0;
            for (std::size_t k = 0; k < C; ++k)
                Uq(a, k) = T(X[k] * pf.D(a, k) / num_traits<T>::from_double(c_i));
        }
    } else if (method == "lin" || method == "gflin" || method == "egflin") {
        // class- or joint-dependent models go to solver_amvald, which handles
        // them; the linearizer kernels do not (solver_amva.m:304-307)
        for (const auto& st : L.stations)
            if (st.cdscaling || st.jdscaling)
                return solver_amvald(L, d, opt, method, converged, init_sol);
        const pfqn::LinearizerMxMethod mx = method == "lin"     ? pfqn::LinearizerMxMethod::Lin
                                            : method == "gflin" ? pfqn::LinearizerMxMethod::Gflin
                                                                : pfqn::LinearizerMxMethod::Egflin;
        if (anyinf) return solver_amvald(L, d, opt, method, converged, init_sol);
        if (allone || opt.multiserver == "krzesinski") {
            std::vector<int> ns;
            for (std::size_t a = 0; a < nq; ++a)
                ns.push_back(std::isfinite(pf.S[a]) ? static_cast<int>(std::llround(pf.S[a])) : 1);
            if (allone) ns.assign(nq, 1);
            const pfqn::LinearizerResult<T> r = pfqn::pfqn_linearizermx(
                pf.lambda, Dm, N, Zm, ns, types, opt.tol, opt.iter_max, mx, Q0);
            Qq = r.Q;
            Uq = r.U;
            X = r.X;
            iters = r.totiter;
        } else if (opt.multiserver == "conway") {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "solver_amva: the Conway multiserver correction evaluates transcendentals, "
                    "which exact arithmetic has no representation for; use the double or real "
                    "backend");
            } else {
                std::vector<int> ns;
                for (std::size_t a = 0; a < nq; ++a)
                    ns.push_back(std::isfinite(pf.S[a]) ? static_cast<int>(std::llround(pf.S[a]))
                                                        : 1);
                const pfqn::LinearizerResult<T> r =
                    pfqn::pfqn_conwayms(pf.D, N, pf.Z, ns, types, opt.tol, opt.iter_max, Q0);
                Qq = r.Q;
                Uq = r.U;
                X = r.X;
                iters = r.totiter;
            }
        } else {
            // default, seidmann, softmin and suri all go to the general AMVA
            return solver_amvald(L, d, opt, method, converged, init_sol);
        }
    } else {
        // multiserver-rule reset rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
        MvaOptions o2 = opt;
        if (o2.multiserver == "conway" || o2.multiserver == "erlang" ||
            o2.multiserver == "krzesinski")
            o2.multiserver = "default";
        return solver_amvald(L, d, o2, method, converged, init_sol);
    }

    Matrix<T> Q(M, C, zero), U(M, C, zero), Tp(M, C, zero), R(M, C, zero);
    for (std::size_t a = 0; a < nq; ++a)
        for (std::size_t c = 0; c < C; ++c) {
            Q(pf.queue_stations[a] - 1, c) = Qq(a, c);
            U(pf.queue_stations[a] - 1, c) = Uq(a, c);
        }
    // Delay stations: Q = X Z with the ORIGINAL think times, for every method.
    // Seidmann folds L(m-1)/m of each multiserver station into Z, but that
    // population is in service at the station and is given back to it just
    // below; charging Zm here too counted it twice and sum(Q) exceeded N.
    for (std::size_t a = 0; a < nz; ++a)
        for (std::size_t c = 0; c < C; ++c) {
            Q(pf.delay_stations[a] - 1, c) = T(X[c] * pf.Z(a, c));
            U(pf.delay_stations[a] - 1, c) = Q(pf.delay_stations[a] - 1, c);
        }
    (void)have_delay_rows;  // the reference overwrites the solved delay rows too
    // Un-apply Seidmann back onto the originating queue. ab and schmidt never
    // had it applied, so they are exempt, which is what the reference tests.
    if (seidmann && !direct_ms) {
        for (std::size_t a = 0; a < nq; ++a) {
            const double c = pf.S[a];
            if (!std::isfinite(c) || c <= 1.0) continue;
            for (std::size_t k = 0; k < C; ++k)
                Q(pf.queue_stations[a] - 1, k) =
                    T(Q(pf.queue_stations[a] - 1, k) +
                      pf.D(a, k) * num_traits<T>::from_double((c - 1.0) / c) * X[k]);
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c) {
            Tp(i, c) = T(d.Vchain(i, c) * X[c]);
            if (Tp(i, c) != zero) R(i, c) = T(Q(i, c) / Tp(i, c));
        }
    // Cycle time excludes the think time actually spent at the delays, which is
    // the ORIGINAL Z summed over them; the Seidmann Zm disagreed with the station
    // residence times sum(R.*V) that Q now reports.
    std::vector<T> Cyc(C, zero);
    for (std::size_t c = 0; c < C; ++c) {
        if (!(X[c] > zero)) continue;
        T z = zero;
        for (std::size_t a = 0; a < nz; ++a) z += pf.Z(a, c);
        Cyc[c] = T(num_traits<T>::from_double(d.Nchain[c]) / X[c] - z);
    }

    MvaSolution<T> out;
    out.method = method;
    out.iter = iters;
    if (L.has_class_switching()) {
        // de-aggregation input rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
        const ClassResults<T> cr =
            sn_deaggregate_chain_results(L, d, Matrix<T>(), Matrix<T>(), R, Tp, X);
        out.Q = cr.Q;
        out.U = cr.U;
        out.R = cr.R;
        out.Tp = cr.Tp;
        out.C = cr.C;
        out.X = cr.X;
    } else {
        out.Q = Q;
        out.U = U;
        out.R = R;
        out.Tp = Tp;
        out.X = X;
        out.C = Cyc;
    }
    return out;
}

// ---------------------------------------------------------------------------
// solver_mvald: the exact load-dependent recursion
// ---------------------------------------------------------------------------

/**
 * Port of `solver_mvald.m`: exact MVA on a load-dependent model, through
 * `pfqn_mvaldmx`.
 *
 * The rate lattice `mu` is built per station: an infinite server gets the
 * linear ramp 1, 2, ..., Nt (which is what makes it a delay in a
 * load-dependent recursion), a station with `lldscaling` gets its lattice, and
 * everything else gets ones.
 *
 * THE UTILIZATION IS NOT the recursion's. Under load-dependent scaling the
 * reference replaces it with the carried load over the EFFECTIVE capacity
 * `max(nservers, max(lldscaling))`, the NC convention, rather than the
 * P(busy)-style estimator `pfqn_mvaldmx` returns; the two disagree wherever the
 * lattice exceeds the server count.
 */
template <class T>
MvaSolution<T> solver_mvald(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const ChainDemands<T> d = sn_get_demands_chain(L);
    const std::size_t M = L.nstations, C = L.nchains;

    double Nt = 0.0;
    for (std::size_t c = 0; c < C; ++c)
        if (std::isfinite(d.Nchain[c])) Nt += d.Nchain[c];
    const std::size_t NT = static_cast<std::size_t>(std::llround(Nt));

    // An open chain contributes its arrival rate and is marked OPEN_CLASS. The
    // rate is read at CHAIN level: STchain at an open chain's reference station
    // (its Source) is one over the SUM of the class arrival rates, which also
    // covers a chain whose classes arrive at several rates.
    std::vector<T> lambda(C, zero);
    std::vector<int> N(C, 0);
    std::vector<bool> openChain(C, false);
    std::size_t nOpenChains = 0;
    for (std::size_t c = 0; c < C; ++c) {
        if (std::isfinite(d.Nchain[c])) {
            N[c] = static_cast<int>(std::llround(d.Nchain[c]));
            continue;
        }
        N[c] = pfqn::OPEN_CLASS;
        openChain[c] = true;
        ++nOpenChains;
        const T st = d.STchain(d.refstatchain[c] - 1, c);
        if (st > zero) lambda[c] = T(one / st);
    }

    Matrix<T> Qchain(M, C, zero), Uchain(M, C, zero);
    std::vector<T> Xchain(C, zero);
    if (nOpenChains == 0) {
        // PURELY CLOSED. Every station enters the recursion, an infinite server
        // as the load-dependent rate mu(n)=n, exact because n cannot then exceed
        // the closed population.
        if (NT == 0) throw UnsupportedError("solver_mvald: the model has no closed population");
        Matrix<T> mu(M, NT, one);
        for (std::size_t i = 0; i < M; ++i) {
            if (std::isinf(L.stations[i].nservers)) {
                for (std::size_t n = 0; n < NT; ++n)
                    mu(i, n) = num_traits<T>::from_int(static_cast<long>(n) + 1);
            } else if (!L.stations[i].lldscaling.empty()) {
                const std::vector<T>& a = L.stations[i].lldscaling;
                for (std::size_t n = 0; n < NT; ++n) mu(i, n) = n < a.size() ? a[n] : a.back();
            }
        }
        Matrix<T> Dm(M, C, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < C; ++c) Dm(i, c) = d.Lchain(i, c);
        const pfqn::MvaResult<T> pf = pfqn::pfqn_mvaldmx(lambda, Dm, N, Matrix<T>(), mu);
        Xchain = pf.XN;
        Qchain = pf.QN;
        Uchain = pf.UN;
    } else {
        // MIXED OR PURELY OPEN. Three kinds of row are not the same thing to
        // pfqn_mvaldmx and have to be separated before it is called; this is the
        // partition solver_ncld makes for the same recursion.
        //  - THE SOURCE IS NOT A STATION. Its chain demand is the interarrival
        //    time 1/lambda, so it carries offered load Lo=1 exactly and
        //    pfqn_ldmx_ec then forms 1/(1-Lo/mu) = inf, poisoning every chain.
        //  - A DELAY IS AN INFINITE SERVER FOR THE OPEN CHAINS TOO. mu(n)=n cut
        //    at the closed population declares it saturated at NT jobs. It enters
        //    as chain think time instead and its queue length is X*L, exact.
        //  - A QUEUEING STATION KEEPS ITS WHOLE RATE ROW. pfqn_ldmx_ec reads the
        //    limited-load-dependence level b off the row itself, so a row cut at
        //    the closed population is read as a slower station, and with no closed
        //    class at all it collapses to mu(1), a single fixed-rate server.
        std::vector<bool> sourceStation(M, false), delayStation(M, false);
        for (std::size_t c = 0; c < C; ++c)
            if (openChain[c]) sourceStation[static_cast<std::size_t>(d.refstatchain[c] - 1)] = true;
        std::vector<std::size_t> queueStations;
        for (std::size_t i = 0; i < M; ++i) {
            if (sourceStation[i]) continue;
            if (std::isinf(L.stations[i].nservers))
                delayStation[i] = true;
            else
                queueStations.push_back(i);
        }
        const std::size_t nq = queueStations.size();

        Matrix<T> Z(1, C, zero);
        for (std::size_t c = 0; c < C; ++c) {
            T z = zero;
            for (std::size_t i = 0; i < M; ++i)
                if (delayStation[i]) z = T(z + d.Lchain(i, c));
            Z(0, c) = z;
        }

        std::size_t ncol = NT > 0 ? NT : 1;
        for (std::size_t k = 0; k < nq; ++k) {
            // first column of the trailing constant run, the level b of pfqn_ldmx_ec
            const std::vector<T>& a = L.stations[queueStations[k]].lldscaling;
            std::size_t b = a.size();
            while (b > 1 && a[b - 2] == a[b - 1]) --b;
            ncol = std::max(ncol, b);
        }
        Matrix<T> mu(nq, ncol, one);
        for (std::size_t k = 0; k < nq; ++k) {
            const std::vector<T>& a = L.stations[queueStations[k]].lldscaling;
            if (a.empty()) continue;
            // held at the row's last value past its own end: limited load dependence
            for (std::size_t n = 0; n < ncol; ++n) mu(k, n) = n < a.size() ? a[n] : a.back();
        }
        Matrix<T> Dq(nq, C, zero);
        for (std::size_t k = 0; k < nq; ++k)
            for (std::size_t c = 0; c < C; ++c) Dq(k, c) = d.Lchain(queueStations[k], c);

        const pfqn::MvaResult<T> pf = pfqn::pfqn_mvaldmx(lambda, Dq, N, Z, mu);
        Xchain = pf.XN;
        for (std::size_t k = 0; k < nq; ++k)
            for (std::size_t c = 0; c < C; ++c) {
                Qchain(queueStations[k], c) = pf.QN(k, c);
                Uchain(queueStations[k], c) = pf.UN(k, c);
            }
        for (std::size_t i = 0; i < M; ++i)
            if (delayStation[i])
                for (std::size_t c = 0; c < C; ++c)
                    // infinite server: X*L for a closed chain, lambda*L for an open one
                    Qchain(i, c) = T(d.Lchain(i, c) * Xchain[c]);
    }

    // Rchain IS THE PER-VISIT RESPONSE TIME, Qchain / Tchain, not the residence
    // time Qchain / Xchain. `sn_deaggregate_chain_results` rebuilds the class
    // queue length as `Rchain * Xchain * Vchain(i,c)/Vchain(refstat,c)`, so it
    // multiplies the visit ratio back IN; handing it a residence time counts
    // that ratio twice and the reported queue lengths then no longer sum to the
    // population. Invisible on every load-dependent model whose LD station has
    // the reference station's chain visits (the shipped examples all do), and
    // decisive on one that does not -- the closed delayed-hit retrieval cache,
    // whose fetch station is visited once per MISS. Every sibling analyzer here
    // (:300, :1046, :1690) already divides by Tchain; this line did not.
    // MATLAB's solver_mvald.m:45 carried the same divisor and is fixed with it;
    // Python's exact LD path already divided by Tchain and was right.
    Matrix<T> Tchain(M, C, zero), Rchain(M, C, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c) {
            Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));
            if (Tchain(i, c) > zero) Rchain(i, c) = T(Qchain(i, c) / Tchain(i, c));
        }

    const ClassResults<T> cr =
        sn_deaggregate_chain_results(L, d, Matrix<T>(), Uchain, Rchain, Tchain, Xchain);
    MvaSolution<T> out;
    out.Q = cr.Q;
    out.U = cr.U;
    out.R = cr.R;
    out.Tp = cr.Tp;
    out.C = cr.C;
    out.X = cr.X;
    out.method = "exact";
    out.iter = 1;

    for (std::size_t i = 0; i < M; ++i) {
        if (L.stations[i].lldscaling.empty() || !std::isfinite(L.stations[i].nservers)) continue;
        double ceff = L.stations[i].nservers;
        for (const T& v : L.stations[i].lldscaling)
            ceff = std::max(ceff, num_traits<T>::to_double(v));
        for (std::size_t r = 0; r < L.nclasses; ++r) {
            if (L.disabled[i][r]) continue;
            const T st = L.rates(i, r) > zero ? T(one / L.rates(i, r)) : zero;
            if (st > zero) out.U(i, r) = T(out.Tp(i, r) * st / num_traits<T>::from_double(ceff));
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// solver_mva_sum: the summation method
// ---------------------------------------------------------------------------

/**
 * Port of `solver_mva_sum.m`: the SUM / ESUM summation method (Bolch et al.,
 * Secs. 9.2, 10.1.4.4 and 10.1.5).
 *
 * A closed model goes to `sum_closed`; an open or mixed one to `sum_closing`
 * with the reference's Kclosed = 5000. The SCV each station is given is the
 * ESUM discrimination: FCFS and SIRO are service-time sensitive and get the
 * chain SCV, while PS, LCFSPR and the infinite servers are insensitive and are
 * passed 1 -- handing them their real SCV would apply a correction that the
 * product form says does not exist.
 */
template <class T>
MvaSolution<T> solver_mva_sum(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using qn::SchedStrategy;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const ChainDemands<T> d = sn_get_demands_chain(L);
    const std::size_t M = L.nstations, K = L.nchains;

    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mva_sum: the summation method locates the root of its population constraint "
            "by bisection and evaluates the Erlang-C waiting probability, neither of which exact "
            "arithmetic has a representation for; use the double or real backend");
    } else {
    std::vector<std::size_t> rows;
    std::vector<sum::Servers> mi;
    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy sc = L.stations[i].sched;
        if (sc == SchedStrategy::EXT) continue;  // the external world is lambda
        if (sc == SchedStrategy::INF) {
            rows.push_back(i);
            mi.push_back(sum::Servers::inf());
            continue;
        }
        if (sc == SchedStrategy::PS || sc == SchedStrategy::LCFSPR ||
            sc == SchedStrategy::FCFS || sc == SchedStrategy::SIRO) {
            rows.push_back(i);
            const double c = L.stations[i].nservers;
            mi.push_back(std::isfinite(c) ? sum::Servers::of(static_cast<long>(std::llround(c)))
                                          : sum::Servers::inf());
            continue;
        }
        throw UnsupportedError(std::string("solver_mva_sum: the summation method does not "
                                           "support ") +
                               lang::sched_to_text(sc) + " scheduling");
    }

    const std::size_t Mr = rows.size();
    Matrix<T> Lm(Mr, K, zero), scv(Mr, K, one);
    for (std::size_t a = 0; a < Mr; ++a) {
        const std::size_t i = rows[a];
        const bool sensitive = L.stations[i].sched == SchedStrategy::FCFS ||
                               L.stations[i].sched == SchedStrategy::SIRO;
        for (std::size_t c = 0; c < K; ++c) {
            Lm(a, c) = T(d.STchain(i, c) * d.Vchain(i, c));
            if (!sensitive) continue;
            const double v = num_traits<T>::to_double(d.SCVchain(i, c));
            if (std::isfinite(v) && v > 0.0) scv(a, c) = num_traits<T>::from_double(v);
        }
    }

    std::vector<long> N(K, 0);
    std::vector<T> Z(K, zero);
    std::vector<std::size_t> ocl;
    for (std::size_t c = 0; c < K; ++c) {
        if (std::isfinite(d.Nchain[c])) N[c] = static_cast<long>(std::llround(d.Nchain[c]));
        else ocl.push_back(c);
    }

    Matrix<T> Qrows, Urows;
    std::vector<T> Xchain;
    std::size_t iters = 0;
    if (ocl.empty()) {
        sum::SumOptions so;
        so.tol = opt.iter_tol;
        so.maxiter = static_cast<std::size_t>(opt.iter_max);
        const sum::SumClosedResult<T> r = sum::sum_closed(Lm, N, Z, mi, scv, so);
        Qrows = r.QN;
        Urows = r.UN;
        Xchain = r.XN;
        iters = r.it;
    } else {
        std::vector<T> lambda(K, zero), scva(K, one);
        for (std::size_t c : ocl) {
            const T st = d.STchain(d.refstatchain[c] - 1, c);
            // AN OPEN CHAIN WITH NO ARRIVAL RATE IS NOT A CHAIN OF RATE ZERO.
            // Leaving lambda[c] at zero makes `sum_closing` treat it as a CLOSED
            // chain of zero population, which is a different model solved
            // without complaint. The service time at the reference station is
            // where the arrival rate comes from, so a non-positive one means the
            // chain is unspecified, not idle.
            if (!(st > zero))
                throw InputError(
                    "solver_mva_sum: open chain " + std::to_string(c + 1) +
                    " has a non-positive service time at its reference station, so it carries no "
                    "arrival rate; the summation method cannot place it");
            lambda[c] = T(one / st);
            const double v = num_traits<T>::to_double(d.SCVchain(d.refstatchain[c] - 1, c));
            if (std::isfinite(v) && v > 0.0) scva[c] = num_traits<T>::from_double(v);
        }
        sum::ClosingOptions co;
        co.sum.tol = opt.iter_tol;
        co.sum.maxiter = static_cast<std::size_t>(opt.iter_max);
        const sum::SumClosingResult<T> r = sum::sum_closing(lambda, scva, Lm, mi, scv, N, Z, co);
        Qrows = r.QN;
        Urows = r.UN;
        Xchain = r.XN;
        iters = r.it;
    }

    Matrix<T> Qchain(M, K, zero), Uchain(M, K, zero), Tchain(M, K, zero), Rchain(M, K, zero);
    for (std::size_t a = 0; a < Mr; ++a)
        for (std::size_t c = 0; c < K; ++c) {
            Qchain(rows[a], c) = Qrows(a, c);
            Uchain(rows[a], c) = Urows(a, c);
        }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));
            if (Tchain(i, c) > zero) Rchain(i, c) = T(Qchain(i, c) / Tchain(i, c));
        }
    for (std::size_t c = 0; c < K; ++c) {
        if (d.Nchain[c] != 0.0) continue;
        Xchain[c] = zero;
        for (std::size_t i = 0; i < M; ++i) {
            Qchain(i, c) = zero;
            Uchain(i, c) = zero;
            Rchain(i, c) = zero;
            Tchain(i, c) = zero;
        }
    }

    const ClassResults<T> cr =
        sn_deaggregate_chain_results(L, d, Matrix<T>(), Matrix<T>(), Rchain, Tchain, Xchain);
    MvaSolution<T> out;
    out.Q = cr.Q;
    out.U = cr.U;
    out.R = cr.R;
    out.Tp = cr.Tp;
    out.C = cr.C;
    out.X = cr.X;
    out.method = opt.method;
    out.iter = static_cast<int>(iters);
    return out;
    }
}

// ---------------------------------------------------------------------------
// Blocking-after-service: the SQD handler
// ---------------------------------------------------------------------------

/**
 * Port of `solver_sqd.m`: Smith queue decomposition for blocking-after-service.
 *
 * `npfqn_sqd` is single-chain by construction -- it solves ONE circulating
 * population -- so a multichain model is refused. The reference warns and
 * returns a matrix of NaN there; refusing by name is the same information
 * without a result that reads as solved.
 *
 * The unpacking is the one `npfqn_sqd.h` documents: the chain-aggregated demand
 * and visit columns, an infinite capacity at every delay, and `sn.rt` with its
 * stateful indexing.
 */
template <class T>
MvaSolution<T> solver_sqd(const qn::NetworkStruct<T>& L, const ChainDemands<T>& d) {
    using qn::SchedStrategy;
    // The effective-rate calibration evaluates exp, log and real powers, and the
    // M/M/1/K blocking probability a real power, so there is no exact-field
    // version of this method. Refuse at RUN time rather than at compile time,
    // as every other transcendental analyzer in this port does.
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)d;
        throw UnsupportedError(
            "solver_sqd: the blocking-after-service calibration evaluates exp, log and real "
            "powers, which exact rational arithmetic cannot represent");
    } else {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations;
    if (L.nchains != 1)
        throw UnsupportedError(
            "solver_sqd: Smith queue decomposition supports single-chain closed networks only; "
            "this model has " + std::to_string(L.nchains) + " chains");
    // AND CLOSED, which the chain count alone does not establish. The automatic
    // route here is gated by `mva_is_bas_model`, which excludes open classes, but
    // the explicit `method == "sqd"` route is not: an open model reaches
    // npfqn_sqd with a population counting only the closed classes and gets a
    // silently wrong answer rather than a refusal.
    if (L.has_open_classes())
        throw UnsupportedError(
            "solver_sqd: Smith queue decomposition is defined for a CLOSED population; this model "
            "has an open class, whose jobs the fixed point cannot count");

    std::vector<T> ST(M, zero), V(M, zero), cap(M, zero);
    std::vector<bool> isDelay(M, false);
    std::vector<std::size_t> s2sf(M, 0);
    const T inf = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < M; ++i) {
        ST[i] = d.STchain(i, 0);
        V[i] = d.Vchain(i, 0);
        isDelay[i] = (L.stations[i].sched == SchedStrategy::INF ||
                      L.stations[i].sched == SchedStrategy::EXT);
        // a capacity above 1e14 is MATLAB's own "effectively unbounded" test
        const double c = L.cap.empty() ? std::numeric_limits<double>::infinity() : L.cap[i];
        cap[i] = (isDelay[i] || !(c <= 1e14)) ? inf : num_traits<T>::from_double(c);
        s2sf[i] = L.stateful_of_station(i + 1);
    }

    npfqn::SqdOptions<T> sopt;
    const npfqn::SqdResult<T> r =
        npfqn::npfqn_sqd(ST, V, cap, isDelay, L.rt, s2sf, L.nclasses,
                         static_cast<int>(std::llround(L.nclosedjobs())), sopt);

    // Vchain is normalized to 1 at the reference station, so the per-station
    // throughput there IS the chain reference throughput.
    const std::size_t refstat = d.refstatchain.empty() ? 1 : d.refstatchain[0];
    Matrix<T> Tchain(M, 1, zero), Qchain(M, 1, zero), Uchain(M, 1, zero), Rchain(M, 1, zero);
    for (std::size_t i = 0; i < M; ++i) {
        Tchain(i, 0) = r.X[i];
        Qchain(i, 0) = r.Q[i];
        Uchain(i, 0) = r.U[i];
        Rchain(i, 0) = r.R[i];
    }
    std::vector<T> Xchain(1, r.X[refstat - 1]);

    const ClassResults<T> cr =
        sn_deaggregate_chain_results(L, d, Qchain, Uchain, Rchain, Tchain, Xchain);
    MvaSolution<T> out;
    out.Q = cr.Q;
    out.U = cr.U;
    out.R = cr.R;
    out.Tp = cr.Tp;
    out.C = cr.C;
    out.X = cr.X;
    out.method = "sqd";
    out.lG = std::numeric_limits<double>::quiet_NaN();
    return out;
    }
}

/**
 * Port of `isBasModel` in `solver_mva_analyzer.m`: a closed single-CHAIN model
 * with at least one blocking-after-service drop rule, which SQD handles and
 * neither exact MVA nor AMVA does. This is the DISPATCH test.
 */
template <class T>
bool mva_is_bas_model(const qn::NetworkStruct<T>& L) {
    if (L.nchains != 1 || !(L.nclosedjobs() > 0.0) || L.droprule.empty()) return false;
    for (const auto& c : L.classes)
        if (std::isinf(c.population)) return false;  // an open class present
    for (const auto& row : L.droprule)
        for (lang::DropStrategy s : row)
            if (s == lang::DropStrategy::BAS) return true;
    return false;
}

/**
 * Port of `matlab/src/api/sn/sn_is_bas_model.m`: a closed single-CLASS model
 * with a BAS drop rule. This is the GATE test, and it is deliberately NARROWER
 * than mva_is_bas_model above: the reference keeps two predicates because SQD's
 * Smith decomposition models ONE circulating population, so a chain built from
 * several classes is dispatched to SQD but is NOT exempted from the
 * finite-capacity refusal. Collapsing the two would let a multiclass BAS model
 * past the gate on the strength of a decomposition that does not represent it.
 */
template <class T>
bool sn_is_bas_model(const qn::NetworkStruct<T>& L) {
    if (L.nclasses != 1 || !(L.nclosedjobs() > 0.0) || L.droprule.empty()) return false;
    for (const auto& c : L.classes)
        if (std::isinf(c.population)) return false;  // an open class present
    for (const auto& row : L.droprule)
        for (lang::DropStrategy s : row)
            if (s == lang::DropStrategy::BAS) return true;
    return false;
}

/**
 * Port of `matlab/src/api/sn/sn_is_mm1k_loss.m`: a single-class open
 * Source-Queue-Sink system whose queue is a single-server exponential M/M/1/K
 * with tail drop. Both the MVA moment-based branch (qsys_mg1k_loss_mgs) and the
 * NC probability-based one (qsys_mm1k_loss) gate their finite-capacity loss
 * path on it, and it is the one shape whose finite buffer IS honoured here.
 */
template <class T>
bool sn_is_mm1k_loss(const qn::NetworkStruct<T>& L) {
    if (L.nclasses != 1 || L.nclosedjobs() != 0.0 || L.nodes.size() != 3) return false;
    std::size_t nsrc = 0, nq = 0, nsnk = 0;
    for (const qn::NodeDef& nd : L.nodes) {
        if (nd.nodetype == qn::NodeType::Source) ++nsrc;
        else if (nd.nodetype == qn::NodeType::Queue) ++nq;
        else if (nd.nodetype == qn::NodeType::Sink) ++nsnk;
    }
    if (nsrc != 1 || nq != 1 || nsnk != 1) return false;
    std::size_t qist = 0, sist = 0;
    for (std::size_t i = 0; i < L.stations.size(); ++i) {
        if (L.stations[i].nodetype == qn::NodeType::Queue) qist = i + 1;
        else if (L.stations[i].nodetype == qn::NodeType::Source) sist = i + 1;
    }
    if (qist == 0 || sist == 0) return false;
    if (L.stations[qist - 1].nservers != 1.0) return false;
    if (L.droprule.size() < qist || L.droprule[qist - 1].empty()) return false;
    if (L.droprule[qist - 1][0] != lang::DropStrategy::DROP) return false;
    if (L.cap.size() < qist) return false;
    if (!std::isfinite(L.cap[qist - 1]) || L.cap[qist - 1] <= 0.0) return false;
    const double scvs = num_traits<T>::to_double(L.scv(sist - 1, 0));
    const double scvq = num_traits<T>::to_double(L.scv(qist - 1, 0));
    return std::fabs(scvs - 1.0) <= 1e-6 && std::fabs(scvq - 1.0) <= 1e-6;
}

/**
 * Conservative form of the branch test in solver_amva: true when this model may be solved by
 * the product-form AMVA kernels rather than by solver_amvald. The mixed case is reported true
 * for every resolved method, while the branch itself takes it only for some, so a caller is
 * never told that solver_amvald will run when it might not.
 */
template <class T>
bool amva_uses_pf_kernels(const qn::NetworkStruct<T>& L) {
    bool has_ld = false;
    for (const auto& st : L.stations)
        if (!st.lldscaling.empty() || st.cdscaling || st.jdscaling) has_ld = true;
    return L.has_product_form_not_het_fcfs() && !has_ld &&
           (!L.has_open_classes() || L.has_product_form());
}

/**
 * True when the MVA path this model already dispatches to carries a class-level interlock
 * matrix (Franks 1999, Eq. 4.7) itself, so that supplying one does not silently move the model
 * to a DIFFERENT algorithm.
 *
 * Only two kernels implement the correction: pfqn_mva (exact, closed single-server) and the
 * AMVA forward step of solver_amvald. A model that would otherwise be solved by exact
 * multiserver or mixed MVA, or by the product-form AMVA kernels, cannot take the matrix
 * without swapping its algorithm, and the swap is worth far more than the correction it
 * carries: inside SolverLN it can turn a converging Picard iteration into a limit cycle. A
 * caller holding a matrix such a model cannot carry must apply its own correction instead.
 */
template <class T>
bool mva_carries_interlock(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    const std::string& method = opt.method;
    bool has_open = false, has_closed = false, integral = true, has_finite_server = false;
    double maxsrv = -1.0, Nsum = 0.0;
    for (const auto& c : L.classes) {
        if (std::isinf(c.population)) {
            has_open = true;
        } else {
            if (c.population > 0.0) has_closed = true;
            if (c.population != std::floor(c.population)) integral = false;
        }
        Nsum += c.population;
    }
    for (const auto& st : L.stations)
        if (std::isfinite(st.nservers)) {
            has_finite_server = true;
            maxsrv = std::max(maxsrv, st.nservers);
        }
    // pfqn_mva takes the matrix for a closed single-server model, and for nothing else
    const bool pfqn_mva_can_take_it = !has_open && integral && (maxsrv < 0.0 || maxsrv <= 1.0);

    if (method == "exact" || method == "mva") return pfqn_mva_can_take_it;
    static const char* amva_methods[] = {
        "amva",  "bs",     "qd",   "qli",  "fli",   "lin",          "qdlin", "sqni",
        "gflin", "egflin", "ab",   "schmidt", "schmidt-ext", "tay", "scat",  "aql",
        "qsa",   "lcp",    "chow", "pamb", "pami",  "pamt",         "clust", "dmlin",
        "priomva"};
    for (const char* m : amva_methods)
        if (method == m) return !amva_uses_pf_kernels(L);  // only solver_amvald carries it
    if (method == "default") {
        if (mva_is_bas_model(L)) return false;  // solver_sqd has no interlock term
        const bool exact_mixed = has_open && has_closed && has_finite_server && maxsrv == 1.0 &&
                                 L.has_product_form() && integral;
        const bool exact_small = L.nchains <= 4 && Nsum <= 20.0 && L.has_product_form() &&
                                 !L.has_fractional_populations();
        if (exact_mixed || exact_small) return pfqn_mva_can_take_it;
        return !amva_uses_pf_kernels(L);
    }
    // mvac, sqd, sum, qna, rqna and rqt reach neither kernel
    return false;
}

// ---------------------------------------------------------------------------
// solver_mva_analyzer: the `default` ladder
// ---------------------------------------------------------------------------

/** Port of solver_mva_analyzer.m, `default` and the explicit method names. */
template <class T>
MvaSolution<T> solver_mva_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt,
                                   const Matrix<T>& init_sol) {
    const ChainDemands<T> d = sn_get_demands_chain(L);
    bool converged = true;

    if (opt.method == "exact" || opt.method == "mva") return solver_mva(L, d, opt);
    if (opt.method == "sum" || opt.method == "esum") return solver_mva_sum(L, opt);
    if (opt.method == "sqd") return solver_sqd(L, d);
    // `qna` is handled one level up, in mva_dispatch, purely so that
    // solver_qna.h need not include this header; reaching it here means the
    // analyzer was called directly, which no ported path does.
    if (opt.method == "qna")
        throw UnsupportedError(
            "solver_mva_analyzer: 'qna' is dispatched by mva_dispatch, not by the analyzer");
    if (opt.method == "rqt")
        throw UnsupportedError(
            "solver_mva_analyzer: 'rqt' is dispatched by mva_dispatch, not by the analyzer");
    if (opt.method == "rqna")
        throw UnsupportedError(
            "solver_mva_analyzer: 'rqna' is dispatched by mva_dispatch, not by the analyzer");
    if (opt.method != "default") {
        MvaOptions o = opt;
        return solver_amva(L, d, o, init_sol, converged);
    }

    // BAS: neither exact MVA nor AMVA represents the blocked-server time, so the
    // reference sends it to SQD before every product-form test below.
    if (mva_is_bas_model(L)) return solver_sqd(L, d);

    // Force AMVA for class- or joint-dependent models: exact MVA does not
    // support them (solver_mva_analyzer.m:63-67).
    for (const auto& st : L.stations)
        if (st.cdscaling || st.jdscaling) {
            MvaOptions o = opt;
            return solver_amva(L, d, o, init_sol, converged);
        }

    // An interlock matrix that exact MVA cannot honour sends the model to AMVA, which
    // applies the same Eq. (4.7) correction to the arrival-instant queue length: pfqn_mva
    // carries it for closed single-server models only.
    bool il_needs_amva = false;
    if (!opt.interlock.empty()) {
        il_needs_amva = L.has_open_classes();
        for (const auto& st : L.stations)
            if (std::isfinite(st.nservers) && st.nservers > 1.0) il_needs_amva = true;
    }

    // mixed open+closed exact-MVA gate rationale: see _kb/06-solver-catalog.md (cpp port notes: solver_mva.h)
    if (!il_needs_amva && L.has_open_classes()) {
        bool anyopen = false, anyclosed = false, integral = true;
        for (const auto& cl : L.classes) {
            if (std::isinf(cl.population)) {
                anyopen = true;
                continue;
            }
            if (cl.population > 0.0) anyclosed = true;
            if (cl.population != std::floor(cl.population)) integral = false;
        }
        bool anyfinite = false, allone = true;
        for (const auto& st : L.stations)
            if (std::isfinite(st.nservers)) {
                anyfinite = true;
                if (st.nservers != 1.0) allone = false;
            }
        if (anyopen && anyclosed && anyfinite && allone && integral && L.has_product_form())
            return solver_mva(L, d, opt);
    }

    double Nsum = 0.0;
    for (std::size_t c = 0; c < L.nchains; ++c) Nsum += d.Nchain[c];
    if (!il_needs_amva && L.nchains <= 4 && Nsum <= 20.0 && L.has_product_form() &&
        !L.has_fractional_populations())
        return solver_mva(L, d, opt);

    MvaOptions o = opt;
    return solver_amva(L, d, o, init_sol, converged);
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_H
