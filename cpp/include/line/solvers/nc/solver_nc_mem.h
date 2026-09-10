/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_MEM_H
#define LINE_SOLVERS_NC_SOLVER_NC_MEM_H

/**
 * Port of `solver_nc_mem.m` and `solver_nc_mem_supports.m`: the Maximum Entropy
 * Method of Kouvatsos (1994).
 *
 * WHY MEM IS IN A NORMALIZING-CONSTANT SOLVER AT ALL. It is not a
 * normalizing-constant algorithm and shares no code with one. It is here
 * because it answers the question the product-form path cannot: a network whose
 * arrival or service processes are NOT exponential. The product-form analyzer
 * silently exponentializes such a model -- it reads only the mean rate -- while
 * MEM carries the second moment through GE (generalised exponential) building
 * blocks and reports what the variability does.
 *
 * THIS FILE IS A DISPATCHER, NOT AN ALGORITHM. The four algorithms live in
 * `api/me/`, and three of them were already ported:
 *
 *   open      Section 3.2, GE/GE/1, GE/GE/c and GE/GE/inf blocks -> me_oqn
 *   closed    Section 3.3, two-stage pseudo-open plus convolution -> me_cqn
 *   mixed     the two composed by product-form-style conditioning -> me_mqn
 *   blocking  Section 4.1, censored GE/GE/c/0;N blocks           -> me_oqn_blk
 *
 * All four are now ported. The blocking one additionally covers TRANSFER
 * BLOCKING (BAS) through the holding-node expansion of Tahilramani, Manjunath
 * and Bose (1999); see `api/me/me_oqn_blk.h` for why that expansion is needed
 * at all.
 *
 * THE GATE IS PART OF THE METHOD. `solver_nc_mem_supports` is not a
 * convenience: MEM's building blocks are defined for particular network shapes,
 * and a model outside them has no MEM answer at all. The reference returns a
 * REASON naming the first violated rule, and that reason is carried into the
 * exception here so a refusal says which rule and which station.
 *
 * ARITHMETIC. The GE blocks are transcendental throughout; a non-transcendental
 * backend is refused by name.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/sn/sn_get_buffer_size.h"
#include "line/api/me/me_cqn.h"
#include "line/api/me/me_mqn.h"
#include "line/api/me/me_oqn.h"
#include "line/api/me/me_oqn_blk.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/nc/solver_nc.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** The verdict of `solver_nc_mem_supports`. */
struct MemSupport {
    bool supported = false;
    std::string reason;    ///< names the first violated rule; empty when supported
    bool blocking = false; ///< the model carries a binding finite station buffer
};

namespace detail {

/**
 * Kendall's K of a station. Delegates to line::sn::sn_get_buffer_size, which is
 * the one implementation; kept here as a name the NC code already uses.
 */
template <class T>
double buffer_size(const qn::NetworkStruct<T>& sn, std::size_t ist) {
    return line::sn::sn_get_buffer_size(sn, ist);
}

/** The discipline name, for the reason strings. */
inline const char* sched_text(qn::SchedStrategy s) {
    switch (s) {
        case qn::SchedStrategy::INF: return "INF";
        case qn::SchedStrategy::FCFS: return "FCFS";
        case qn::SchedStrategy::PS: return "PS";
        case qn::SchedStrategy::SIRO: return "SIRO";
        case qn::SchedStrategy::LCFS: return "LCFS";
        case qn::SchedStrategy::LCFSPR: return "LCFSPR";
        case qn::SchedStrategy::EXT: return "EXT";
        case qn::SchedStrategy::HOL: return "HOL";
        case qn::SchedStrategy::DPS: return "DPS";
        case qn::SchedStrategy::GPS: return "GPS";
        default: return "an unsupported strategy";
    }
}

}  // namespace detail

/**
 * Port of `solver_nc_mem_supports.m`.
 *
 * @param sn the refreshed struct
 * @return whether MEM applies, the first violated rule if not, and whether the
 *         model needs the finite-buffer (blocking) algorithm
 */
template <class T>
MemSupport solver_nc_mem_supports(const qn::NetworkStruct<T>& sn) {
    MemSupport out;
    const std::size_t M = sn.nstations, R = sn.nclasses;

    bool anyOpen = false, anyClosed = false;
    for (const qn::JobClass& c : sn.classes)
        (std::isinf(c.population) ? anyOpen : anyClosed) = true;
    const bool isopen = anyOpen && !anyClosed;
    const bool isclosed = anyClosed && !anyOpen;
    const bool ismixed = anyOpen && anyClosed;

    // Node types: open models are Source/Queue/Delay/Sink, closed are
    // Queue/Delay.
    for (const qn::NodeDef& nd : sn.nodes) {
        switch (nd.nodetype) {
            case qn::NodeType::Queue:
            case qn::NodeType::Delay:
                break;
            case qn::NodeType::Source:
            case qn::NodeType::Sink:
                if (isclosed) {
                    out.reason = "MEM supports only Queue and Delay nodes in closed models.";
                    return out;
                }
                break;
            default:
                out.reason = "MEM supports only Source, Queue, Delay and Sink nodes.";
                return out;
        }
    }

    if (sn.has_class_switching()) {
        out.reason = "MEM does not support class switching.";
        return out;
    }

    // An absorbing self-loop makes the routing reducible and the geometric
    // feedback transform 1/(1-p_ii) degenerate.
    for (std::size_t ist = 1; ist <= M; ++ist) {
        const std::size_t nd = sn.node_of_station(ist);
        for (std::size_t r = 0; r < R; ++r)
            if (num_traits<T>::to_double(sn.rtnodes((nd - 1) * R + r, (nd - 1) * R + r)) >=
                1.0 - 1e-9) {
                out.reason =
                    "MEM does not support absorbing self-loop routing (reducible network).";
                return out;
            }
    }

    // The PR/HOL constraint formulae are not given in the reference.
    for (std::size_t i = 0; i < M; ++i) {
        const qn::SchedStrategy s = sn.stations[i].sched;
        if (s != qn::SchedStrategy::EXT && s != qn::SchedStrategy::INF &&
            s != qn::SchedStrategy::FCFS && s != qn::SchedStrategy::PS &&
            s != qn::SchedStrategy::SIRO && s != qn::SchedStrategy::LCFS &&
            s != qn::SchedStrategy::LCFSPR) {
            out.reason = std::string("MEM does not support the ") + detail::sched_text(s) +
                         " scheduling strategy.";
            return out;
        }
    }

    if (isopen || ismixed) {
        bool hasSource = false;
        for (const qn::NodeDef& nd : sn.nodes)
            if (nd.nodetype == qn::NodeType::Source) hasSource = true;
        if (!hasSource) {
            out.reason = "MEM requires a Source node when open classes are present.";
            return out;
        }
    }
    if (isclosed || ismixed) {
        // Closed classes build on G/G/1 and G/G/inf blocks only.
        for (std::size_t i = 0; i < M; ++i)
            if (std::isfinite(sn.stations[i].nservers) && sn.stations[i].nservers > 1.0) {
                out.reason =
                    "MEM does not support multiserver stations in closed or mixed models.";
                return out;
            }
    }

    // Finite buffers: the censored block is single class, so a binding buffer
    // is admissible only in a single-class open model.
    std::vector<bool> capped(M, false);
    bool anyCapped = false;
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched == qn::SchedStrategy::EXT) continue;
        if (std::isfinite(detail::buffer_size(sn, i + 1))) {
            capped[i] = true;
            anyCapped = true;
        }
    }
    if (anyCapped) {
        if (!isopen) {
            out.reason = "MEM supports finite station buffers only in open models.";
            return out;
        }
        if (R > 1) {
            out.reason =
                "MEM supports finite station buffers only in single-class models: the censored "
                "GE/GE/c/0;N building block is single class.";
            return out;
        }
        for (std::size_t i = 0; i < M; ++i) {
            if (!capped[i]) continue;
            const double c = sn.stations[i].nservers;
            if (!std::isfinite(c) || c < 1.0) {
                out.reason = "MEM cannot apply a finite buffer to the infinite-server station " +
                             std::to_string(i + 1) + ".";
                return out;
            }
            if (sn.stations[i].sched != qn::SchedStrategy::FCFS) {
                out.reason =
                    "MEM supports finite station buffers only under FCFS scheduling; station " +
                    std::to_string(i + 1) + " uses " + detail::sched_text(sn.stations[i].sched) +
                    ".";
                return out;
            }
            const int dr = sn.stations[i].droprule.empty() ? 0 : sn.stations[i].droprule[0];
            if (dr != 0 && dr != static_cast<int>(lang::DropStrategy::DROP) &&
                dr != static_cast<int>(lang::DropStrategy::BAS)) {
                out.reason =
                    "MEM supports the DROP and BAS drop rules at a finite buffer; station " +
                    std::to_string(i + 1) + " uses another rule.";
                return out;
            }
            // The GE distribution is undefined below scv 1, so a hypo-exponential
            // service is refused rather than approximated.
            if (!sn.disabled[i][0] &&
                num_traits<T>::to_double(sn.scv(i, 0)) < 1.0 - 1e-12) {
                out.reason = "MEM with finite buffers needs a service scv of at least 1 at "
                             "station " + std::to_string(i + 1) +
                             ": the GE distribution is not defined below 1.";
                return out;
            }
        }
        for (std::size_t i = 0; i < M; ++i)
            if (sn.stations[i].sched == qn::SchedStrategy::EXT && !sn.disabled[i][0] &&
                num_traits<T>::to_double(sn.scv(i, 0)) < 1.0 - 1e-12) {
                out.reason = "MEM with finite buffers needs an external interarrival scv of at "
                             "least 1: the GE distribution is not defined below 1.";
                return out;
            }
        out.blocking = true;
    }

    out.supported = true;
    return out;
}

/**
 * Port of `solver_nc_mem.m`.
 *
 * @param sn  the refreshed struct
 * @param opt solver controls; `mem_tol`, `mem_maxiter` are read
 */
template <class T>
NcSolution<T> solver_nc_mem(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc_mem: the Maximum Entropy Method is built on GE (generalised exponential) "
            "blocks whose entropy maximisation is transcendental throughout; this backend has "
            "none");
    } else {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        const std::size_t M = sn.nstations, R = sn.nclasses;

        const MemSupport sup = solver_nc_mem_supports(sn);
        if (!sup.supported) throw UnsupportedError("solver_nc_mem: " + sup.reason);

        me::MeOptions mopt;
        mopt.tol = opt.mem_tol;
        mopt.maxiter = opt.mem_maxiter;

        bool anyOpen = false, anyClosed = false;
        for (const qn::JobClass& c : sn.classes)
            (std::isinf(c.population) ? anyOpen : anyClosed) = true;
        const bool isclosed = anyClosed && !anyOpen;
        const bool isopen = anyOpen && !anyClosed;

        // The per-station service rate and scv, zero where a class is not
        // served (the reference reads NaN there and skips it).
        const auto fill_service = [&](const std::vector<std::size_t>& sts, Matrix<T>& mu,
                                      Matrix<T>& Cs, std::vector<long>& c) {
            mu = Matrix<T>(sts.size(), R, zero);
            Cs = Matrix<T>(sts.size(), R, one);
            c.assign(sts.size(), 1);
            for (std::size_t k = 0; k < sts.size(); ++k) {
                const std::size_t i = sts[k] - 1;
                // me_* takes 0 for an infinite server, where MATLAB takes Inf.
                c[k] = std::isinf(sn.stations[i].nservers)
                           ? 0
                           : static_cast<long>(std::llround(sn.stations[i].nservers));
                for (std::size_t r = 0; r < R; ++r) {
                    if (sn.disabled[i][r] || !(sn.rates(i, r) > zero)) continue;
                    mu(k, r) = sn.rates(i, r);
                    if (sn.scv(i, r) > zero) Cs(k, r) = sn.scv(i, r);
                }
            }
        };
        const auto fill_routing = [&](const std::vector<std::size_t>& sts) {
            std::vector<Matrix<T>> P(R, Matrix<T>(sts.size(), sts.size(), zero));
            for (std::size_t r = 0; r < R; ++r)
                for (std::size_t j = 0; j < sts.size(); ++j) {
                    const std::size_t jn = sn.node_of_station(sts[j]);
                    for (std::size_t k = 0; k < sts.size(); ++k) {
                        const std::size_t kn = sn.node_of_station(sts[k]);
                        P[r](j, k) = sn.rtnodes((jn - 1) * R + r, (kn - 1) * R + r);
                    }
                }
            return P;
        };
        const auto fill_insens = [&](const std::vector<std::size_t>& sts) {
            std::vector<char> ins(sts.size(), 0);
            for (std::size_t k = 0; k < sts.size(); ++k) {
                const qn::SchedStrategy s = sn.stations[sts[k] - 1].sched;
                ins[k] = (s == qn::SchedStrategy::PS || s == qn::SchedStrategy::LCFSPR) ? 1 : 0;
            }
            return ins;
        };

        std::vector<long> N(R, 0);
        for (std::size_t r = 0; r < R; ++r)
            N[r] = std::isinf(sn.classes[r].population)
                       ? 0
                       : static_cast<long>(std::llround(sn.classes[r].population));

        out.sol.Q = Matrix<T>(M, R, zero);
        out.sol.U = Matrix<T>(M, R, zero);
        out.sol.R = Matrix<T>(M, R, zero);
        out.sol.Tp = Matrix<T>(M, R, zero);
        out.sol.X.assign(R, zero);
        out.sol.C.assign(R, zero);
        out.actualmethod = "mem";

        if (isclosed) {
            // Every station is a queueing station; there is no Source.
            std::vector<std::size_t> sts(M);
            for (std::size_t i = 0; i < M; ++i) sts[i] = i + 1;
            Matrix<T> mu, Cs;
            std::vector<long> c;
            fill_service(sts, mu, Cs, c);
            std::vector<long> refstat(R, -1);
            for (std::size_t r = 0; r < R; ++r)
                refstat[r] = static_cast<long>(sn.classes[r].refstat) - 1;
            const me::MeResult<T> res =
                me::me_cqn(M, R, N, mu, Cs, fill_routing(sts), c, refstat, fill_insens(sts), mopt);
            out.sol.Q = res.L;
            out.sol.U = res.rho;
            out.sol.R = res.W;
            out.sol.Tp = res.lambda;
            out.sol.X = res.X;
            for (std::size_t r = 0; r < R; ++r)
                if (res.X[r] > zero)
                    out.sol.C[r] = T(num_traits<T>::from_double(
                                         static_cast<double>(N[r])) /
                                     res.X[r]);
            out.sol.iter = static_cast<int>(res.iter);
            return out;
        }

        // Open and mixed both split off the Source and analyze the rest.
        std::size_t sourceIdx = 0;
        for (std::size_t i = 0; i < M; ++i)
            if (sn.stations[i].nodetype == qn::NodeType::Source) sourceIdx = i + 1;
        if (sourceIdx == 0) throw UnsupportedError("solver_nc_mem: MEM requires a Source node");

        std::vector<std::size_t> qs;
        for (std::size_t i = 1; i <= M; ++i)
            if (i != sourceIdx) qs.push_back(i);
        const std::size_t Mq = qs.size();
        Matrix<T> mu, Cs;
        std::vector<long> c;
        fill_service(qs, mu, Cs, c);

        // External arrivals, spread along the Source's routing.
        Matrix<T> lambda0(Mq, R, zero), Ca0(Mq, R, zero);
        const std::size_t sourceNode = sn.node_of_station(sourceIdx);
        for (std::size_t r = 0; r < R; ++r) {
            if (sn.disabled[sourceIdx - 1][r] || !(sn.rates(sourceIdx - 1, r) > zero)) continue;
            if (!std::isinf(sn.classes[r].population)) continue;
            const T extRate = sn.rates(sourceIdx - 1, r);
            T caExt = one;
            if (sn.scv(sourceIdx - 1, r) > zero) caExt = sn.scv(sourceIdx - 1, r);
            for (std::size_t k = 0; k < Mq; ++k) {
                const std::size_t dn = sn.node_of_station(qs[k]);
                const T p = sn.rtnodes((sourceNode - 1) * R + r, (dn - 1) * R + r);
                if (p > zero) {
                    lambda0(k, r) = T(extRate * p);
                    Ca0(k, r) = caExt;
                }
            }
        }

        if (sup.blocking) {
            // Finite buffers: the censored GE/GE/c/0;N blocks, with the
            // holding-node expansion where the drop rule is BAS. Single class
            // by construction -- the gate above rejects a multiclass model with
            // a binding buffer, because the censored block is single class.
            std::vector<long> Nbuf(Mq, 0), c_blk(Mq, 1);
            std::vector<int> blockrule(Mq, 0);
            for (std::size_t k = 0; k < Mq; ++k) {
                const std::size_t i = qs[k] - 1;
                const double nb = detail::buffer_size(sn, qs[k]);
                // 0 marks an unbounded buffer for me_oqn_blk, where MATLAB uses
                // Inf; likewise 0 servers marks an infinite server.
                Nbuf[k] = std::isfinite(nb) ? static_cast<long>(std::llround(nb)) : 0;
                c_blk[k] = std::isinf(sn.stations[i].nservers)
                               ? 0
                               : static_cast<long>(std::llround(sn.stations[i].nservers));
                const int dr = sn.stations[i].droprule.empty() ? 0 : sn.stations[i].droprule[0];
                blockrule[k] = (dr == static_cast<int>(lang::DropStrategy::BAS)) ? 1 : 0;
            }
            std::vector<T> l0(Mq, zero), ca0(Mq, one), mu1(Mq, zero), cs1(Mq, one);
            for (std::size_t k = 0; k < Mq; ++k) {
                l0[k] = lambda0(k, 0);
                ca0[k] = Ca0(k, 0) > zero ? Ca0(k, 0) : one;
                mu1[k] = mu(k, 0);
                cs1[k] = Cs(k, 0);
            }
            Matrix<T> P1(Mq, Mq, zero);
            {
                const std::vector<Matrix<T>> Pall = fill_routing(qs);
                for (std::size_t a = 0; a < Mq; ++a)
                    for (std::size_t b = 0; b < Mq; ++b) P1(a, b) = Pall[0](a, b);
            }
            me::MeBlkOptions bopt;
            bopt.tol = opt.mem_tol;
            bopt.maxiter = opt.mem_maxiter;
            const me::MeBlkResult<T> br =
                me::me_oqn_blk(Mq, l0, ca0, mu1, cs1, P1, c_blk, Nbuf, blockrule, bopt);
            for (std::size_t k = 0; k < Mq; ++k) {
                out.sol.Q(qs[k] - 1, 0) = br.Q[k];
                out.sol.U(qs[k] - 1, 0) = br.U[k];
                out.sol.R(qs[k] - 1, 0) = br.W[k];
                out.sol.Tp(qs[k] - 1, 0) = br.T_[k];
            }
            // The Source emits at its NOMINAL rate; the jobs lost at a full
            // buffer never reach a station, so the carried flow per station is
            // below it by the loss probability at the entry stations. Only a
            // LOSS station discards -- under BAS the job is held, not dropped.
            if (!sn.disabled[sourceIdx - 1][0] && sn.rates(sourceIdx - 1, 0) > zero) {
                const T x = sn.rates(sourceIdx - 1, 0);
                out.sol.Tp(sourceIdx - 1, 0) = x;
                out.sol.X[0] = x;
                T accepted = x;
                for (std::size_t k = 0; k < Mq; ++k)
                    if (l0[k] > zero && blockrule[k] == 0)
                        accepted -= T(l0[k] * br.PBa[k]);
                if (accepted > zero) {
                    T q = zero;
                    for (std::size_t k = 0; k < Mq; ++k) q += out.sol.Q(qs[k] - 1, 0);
                    out.sol.C[0] = T(q / accepted);
                }
            }
            out.sol.iter = static_cast<int>(br.iter);
            out.sol.method = "mem.blocking";
            out.actualmethod = "mem.blocking";
            return out;
        }

        me::MeResult<T> res;
        if (isopen) {
            res = me::me_oqn(Mq, R, lambda0, Ca0, mu, Cs, fill_routing(qs), c, fill_insens(qs),
                             mopt);
        } else {
            std::vector<char> openCls(R, 0);
            std::vector<long> refstat(R, -1);
            for (std::size_t r = 0; r < R; ++r) {
                openCls[r] = std::isinf(sn.classes[r].population) ? 1 : 0;
                if (!openCls[r]) {
                    const auto it = std::find(qs.begin(), qs.end(), sn.classes[r].refstat);
                    if (it != qs.end())
                        refstat[r] = static_cast<long>(it - qs.begin());
                }
            }
            res = me::me_mqn(Mq, R, openCls, lambda0, Ca0, N, mu, Cs, fill_routing(qs), c, refstat,
                             fill_insens(qs), mopt);
        }

        for (std::size_t k = 0; k < Mq; ++k)
            for (std::size_t r = 0; r < R; ++r) {
                out.sol.Q(qs[k] - 1, r) = res.L(k, r);
                out.sol.U(qs[k] - 1, r) = res.rho(k, r);
                out.sol.R(qs[k] - 1, r) = res.W(k, r);
                out.sol.Tp(qs[k] - 1, r) = res.lambda(k, r);
            }

        // An unstable station's utilization is normalized to sum 1, the LINE
        // convention shared with the product-form runner.
        for (std::size_t k = 0; k < Mq; ++k) {
            const std::size_t i = qs[k] - 1;
            if (!std::isfinite(sn.stations[i].nservers)) continue;
            T tot = zero;
            for (std::size_t r = 0; r < R; ++r) tot += out.sol.U(i, r);
            if (num_traits<T>::to_double(tot) > 1.0)
                for (std::size_t r = 0; r < R; ++r) out.sol.U(i, r) = T(out.sol.U(i, r) / tot);
        }

        for (std::size_t r = 0; r < R; ++r) {
            const bool open = std::isinf(sn.classes[r].population);
            if (open) {
                // The Source reports the external rate; the system response
                // time follows by Little's law over the queueing stations.
                const T x = isopen ? (sn.disabled[sourceIdx - 1][r] ? zero
                                                                    : sn.rates(sourceIdx - 1, r))
                                   : res.X[r];
                if (!(x > zero)) continue;
                out.sol.Tp(sourceIdx - 1, r) = x;
                out.sol.X[r] = x;
                T q = zero;
                for (std::size_t k = 0; k < Mq; ++k) q += out.sol.Q(qs[k] - 1, r);
                out.sol.C[r] = T(q / x);
            } else {
                out.sol.X[r] = res.X[r];
                if (res.X[r] > zero)
                    out.sol.C[r] =
                        T(num_traits<T>::from_double(static_cast<double>(N[r])) / res.X[r]);
            }
        }
        out.sol.iter = static_cast<int>(res.iter);
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_MEM_H
