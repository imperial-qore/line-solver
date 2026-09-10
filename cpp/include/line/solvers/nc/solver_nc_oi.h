/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_OI_H
#define LINE_SOLVERS_NC_SOLVER_NC_OI_H

/**
 * Order-independent (OI) and pass-and-swap (P&S) normalizing-constant analysis.
 *
 * Ports `nc_is_oi_model.m`, `nc_is_pas_model.m`, `solver_nc_oi_analyzer.m` and
 * `solver_nc_pas_is_analyzer.m`.
 *
 * An OI station is a load-dependent server whose total rate mu(n) is invariant
 * under permutations of the ordered microstate, so it is a function of the
 * per-class count vector alone. Such a station is product-form under balanced
 * fairness (Bonald and Proutiere 2003) and the whole network's constant is
 * assembled on the population lattice: `pfqn_ncoi` evaluates the OI stations
 * together with the aggregated delay, and every ordinary BCMP station is folded
 * in by lattice convolution of its load-dependent weight table. Mean queue
 * lengths come from the OI functional-server identity of `pfqn_oi_fnc` (Casale,
 * QEST 2006), E[f(n_i)] = G^+/G - 1 with G^+ the convolution of that station's
 * FNC balance function against the full-network table.
 *
 * A P&S station with a NON-EMPTY swap graph is not order-independent: the
 * ordered-state chain is reducible (Comte and Dorsman 2021) and only the
 * recurrent communicating class carries a product form. Its constant has no
 * exact lattice recursion here and is estimated by the auto-normalized
 * importance sampler `pfqn_pas_is`, which reduces to `pfqn_oi_is` when the swap
 * graph is empty.
 *
 * THE TWO ANALYZERS ARE NOT INTERCHANGEABLE and the dispatch order decides
 * which a model gets: a pure-OI tandem on 'default'/'exact' is caught by the
 * exact analyzer first, and only a genuine swap graph (or an explicit 'is' /
 * 'sampling') reaches the sampler.
 *
 * Arithmetic: the OI analyzer needs `log` for lG and `exp`/`lgamma` for the
 * BCMP weight table, so it is guarded on `has_transcendental`. The sampler needs
 * a random stream and is guarded the same way through `pfqn_pas_is`.
 *
 * NOTE ON PRECISION AT `Real<N>`. `oi_ld_table` accumulates its weight in a
 * DOUBLE logarithm, exactly as the reference's `gammaln`/`exp` do, so the BCMP
 * factor is double precision whatever T is. Everything else -- the balanced-
 * fairness fill, the lattice convolution, the FNC identity -- is field
 * arithmetic in T. A higher-precision run therefore gains nothing on the BCMP
 * stations; the table would have to be built as an exact multinomial product to
 * change that, which is a different algorithm from the reference's.
 */

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "line/api/pfqn/pas_swap2order.h"
#include "line/api/pfqn/pfqn_oi_fnc.h"
#include "line/api/pfqn/pfqn_oi_insvc.h"
#include "line/api/pfqn/pfqn_ncoi.h"
#include "line/api/pfqn/pfqn_pas_is.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"

namespace line {
namespace nc {

namespace detail {

/** Column-major lattice descriptor for 0 <= n <= N, as `oi_lattice`. */
inline void oi_lattice(const std::vector<int>& N, std::vector<std::size_t>& shp,
                       std::vector<std::size_t>& stride, std::size_t& total) {
    const std::size_t R = N.size();
    shp.assign(R, 0);
    stride.assign(R, 1);
    total = 1;
    for (std::size_t d = 0; d < R; ++d) shp[d] = static_cast<std::size_t>(N[d]) + 1;
    for (std::size_t d = 1; d < R; ++d) stride[d] = stride[d - 1] * shp[d - 1];
    for (std::size_t d = 0; d < R; ++d) total *= shp[d];
}

/** Decode a 0-based linear lattice index to its count vector, as `oi_sub`. */
inline std::vector<int> oi_sub(std::size_t i, const std::vector<std::size_t>& shp) {
    std::vector<int> n(shp.size(), 0);
    std::size_t li = i;
    for (std::size_t d = 0; d < shp.size(); ++d) {
        n[d] = static_cast<int>(li % shp[d]);
        li /= shp[d];
    }
    return n;
}

/** Flatten a count vector to its 0-based lattice index. */
inline std::size_t oi_idx(const std::vector<int>& n, const std::vector<std::size_t>& stride) {
    std::size_t i = 0;
    for (std::size_t d = 0; d < n.size(); ++d) i += stride[d] * static_cast<std::size_t>(n[d]);
    return i;
}

/**
 * The canonical ordered microstate holding n_r copies of class r, as
 * `oi_microstate`: any ordering evaluates to the same rate at an OI station.
 */
inline std::vector<std::size_t> oi_microstate(const std::vector<int>& n) {
    std::vector<std::size_t> c;
    for (std::size_t r = 0; r < n.size(); ++r)
        for (int a = 0; a < n[r]; ++a) c.push_back(r + 1);
    return c;
}

/** Count vectors kept per memoized rank rate; a miss past it just recomputes. */
const std::size_t OI_RANK_RATE_MEMO_LIMIT = 1u << 17;

/**
 * An OI rank rate on a per-class COUNT vector, MEMOIZED on that vector.
 *
 * Only the importance sampler should use this. It calls the rate 2(ell+1) times
 * per sampled ordering and walks the same prefix occupancies over and over, so
 * rebuilding the microstate -- and re-evaluating the balance function on it --
 * dominated the estimator. The value depends on nothing but the counts, so the
 * memo returns what the rebuild would have returned, bit for bit. The exact
 * enumerating paths visit each count vector ONCE and must keep calling
 * `oi_microstate` directly: there a memo is pure loss.
 */
template <class T>
inline pfqn::OiRateFun<T> oi_rank_rate(
    const std::function<T(const std::vector<std::size_t>&)>& f) {
    std::shared_ptr<std::map<std::vector<int>, T> > memo(new std::map<std::vector<int>, T>());
    return [f, memo](const std::vector<int>& n) -> T {
        const typename std::map<std::vector<int>, T>::const_iterator it = memo->find(n);
        if (it != memo->end()) return it->second;
        const T v = f(oi_microstate(n));
        if (memo->size() < OI_RANK_RATE_MEMO_LIMIT) (*memo)[n] = v;
        return v;
    };
}

/**
 * Forward balanced-fairness fill of an OI balance function, as `oi_phi`:
 * Phi(0) = 1 and Phi(n) = (1/mu(n)) sum_{r : n_r > 0} Phi(n - e_r).
 */
template <class T>
std::vector<T> oi_phi(const std::function<T(const std::vector<int>&)>& oirate,
                      const std::vector<int>& N) {
    std::vector<std::size_t> shp, stride;
    std::size_t total = 0;
    oi_lattice(N, shp, stride, total);
    const std::size_t R = shp.size();
    std::vector<T> Phi(total, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < total; ++i) {
        const std::vector<int> n = oi_sub(i, shp);
        int tot = 0;
        for (int v : n) tot += v;
        if (tot == 0) {
            Phi[i] = num_traits<T>::from_int(1);
            continue;
        }
        T s = num_traits<T>::from_int(0);
        for (std::size_t r = 0; r < R; ++r)
            if (n[r] > 0) s += Phi[i - stride[r]];
        Phi[i] = T(s / oirate(n));
    }
    return Phi;
}

/**
 * The BCMP load-dependent weight table of one station, as `oi_ld_table`:
 * W(n) = |n|!/prod(n_r!) prod_r D_r^{n_r} / prod_{k=1}^{|n|} beta(k), with
 * beta(k) = min(k, c). A class with zero demand contributes nothing, so any
 * lattice point that places a job there has weight zero.
 */
template <class T>
std::vector<T> oi_ld_table(const std::vector<T>& Dq, double c,
                           const std::vector<std::size_t>& shp, std::size_t total) {
    const std::size_t R = shp.size();
    std::vector<T> W(total, num_traits<T>::from_int(0));
    const double cc = (std::isfinite(c) && c > 0.0) ? c : 1.0;
    for (std::size_t i = 0; i < total; ++i) {
        const std::vector<int> n = oi_sub(i, shp);
        int tot = 0;
        for (int v : n) tot += v;
        double logf = std::lgamma(static_cast<double>(tot) + 1.0);
        bool ok = true;
        for (std::size_t r = 0; r < R; ++r)
            if (n[r] > 0) {
                const double d = num_traits<T>::to_double(Dq[r]);
                if (!(d > 0.0)) {
                    ok = false;
                    break;
                }
                logf += n[r] * std::log(d) - std::lgamma(static_cast<double>(n[r]) + 1.0);
            }
        if (!ok) continue;
        for (int k = 1; k <= tot; ++k)
            logf -= std::log(std::min(static_cast<double>(k), cc));
        W[i] = num_traits<T>::from_double(std::exp(logf));
    }
    return W;
}

/** Lattice convolution C(m) = sum_{0<=a<=m} A(a) B(m-a), as `oi_conv`. */
template <class T>
std::vector<T> oi_conv(const std::vector<T>& A, const std::vector<T>& B,
                       const std::vector<std::size_t>& shp,
                       const std::vector<std::size_t>& stride, std::size_t total) {
    const std::size_t R = shp.size();
    std::vector<std::vector<int>> subs(total);
    for (std::size_t i = 0; i < total; ++i) subs[i] = oi_sub(i, shp);
    std::vector<T> C(total, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < total; ++i) {
        T acc = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j <= i; ++j) {
            bool le = true;
            for (std::size_t d = 0; d < R && le; ++d)
                if (subs[j][d] > subs[i][d]) le = false;
            if (!le) continue;
            std::size_t off = 0;
            for (std::size_t d = 0; d < R; ++d)
                off += stride[d] * static_cast<std::size_t>(subs[i][d] - subs[j][d]);
            acc += T(A[j] * B[off]);
        }
        C[i] = acc;
    }
    return C;
}

/** G^+ = sum_{0<=b<=N} Psi(b) G(N-b), as `oi_fnc_mean`. */
template <class T>
T oi_fnc_mean(const std::vector<T>& Psi, const std::vector<T>& Gfull,
              const std::vector<std::size_t>& shp, const std::vector<std::size_t>& stride,
              std::size_t total) {
    const T zero = num_traits<T>::from_int(0);
    T val = zero;
    for (std::size_t i = 0; i < total; ++i) {
        if (Psi[i] == zero) continue;
        const std::vector<int> b = oi_sub(i, shp);
        std::size_t off = 0;
        for (std::size_t d = 0; d < shp.size(); ++d)
            off += stride[d] * (shp[d] - 1 - static_cast<std::size_t>(b[d]));
        val += T(Psi[i] * Gfull[off]);
    }
    return val;
}

/** Per-class visits normalized to the reference station, as both analyzers do. */
template <class T>
Matrix<T> oi_visits(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations, K = sn.nclasses;
    Matrix<T> V(M, K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        std::size_t chain = 0;
        for (std::size_t c = 0; c < sn.nchains; ++c)
            if (sn.chains[c][r]) chain = c + 1;
        if (chain == 0) continue;
        for (std::size_t i = 0; i < M; ++i)
            V(i, r) = sn.visits[chain - 1](sn.stateful_of_station(i + 1) - 1, r);
        const T vref = V(sn.classes[r].refstat - 1, r);
        if (vref > zero)
            for (std::size_t i = 0; i < M; ++i) V(i, r) = T(V(i, r) / vref);
    }
    return V;
}

/** True when the FCFS / SIRO rates at a station vary across populated classes. */
template <class T>
bool class_dependent_fcfs_rate(const qn::NetworkStruct<T>& sn, std::size_t i) {
    double lo = 0.0, hi = 0.0;
    bool any = false;
    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        if (!(sn.classes[r].population > 0.0)) continue;
        const double v = num_traits<T>::to_double(sn.rates(i, r));
        if (!std::isfinite(v)) continue;
        if (!any) {
            lo = hi = v;
            any = true;
        } else {
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
    }
    return any && (hi - lo) > 1e-9 * hi;
}

}  // namespace detail

/**
 * Port of `nc_is_oi_model.m`: a closed network with at least one OI station and
 * nothing but BCMP product-form stations besides.
 *
 * The OI requirement is what keeps a pure-BCMP network on the ordinary (faster)
 * normalizing-constant path rather than the lattice one.
 */
template <class T>
bool nc_is_oi_model(const qn::NetworkStruct<T>& sn) {
    using qn::SchedStrategy;
    for (const qn::JobClass& c : sn.classes)
        if (std::isinf(c.population)) return false;
    bool hasOI = false;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const SchedStrategy s = sn.stations[i].sched;
        if (s == SchedStrategy::INF) continue;
        if (s == SchedStrategy::PAS || s == SchedStrategy::OI) {
            if (!qn::station_swap_graph_is_zero(sn, i + 1)) return false;  // genuine P&S
            hasOI = true;
        } else if (s == SchedStrategy::PS || s == SchedStrategy::LCFSPR ||
                   s == SchedStrategy::SIRO || s == SchedStrategy::FCFS) {
            if ((s == SchedStrategy::FCFS || s == SchedStrategy::SIRO) &&
                detail::class_dependent_fcfs_rate(sn, i))
                return false;
        } else {
            return false;  // an unsupported (non-product-form) station
        }
    }
    return hasOI;
}

/**
 * Port of `nc_is_pas_model.m`: a closed two-station OI / P&S tandem.
 *
 * The swap graph may be empty, since an OI queue is exactly the P&S
 * specialization with a zero graph and `pfqn_pas_is` reduces to `pfqn_oi_is`
 * there; the exact analyzer is what keeps a pure-OI tandem off this path on
 * 'default' and 'exact'.
 */
template <class T>
bool nc_is_pas_model(const qn::NetworkStruct<T>& sn) {
    using qn::SchedStrategy;
    for (const qn::JobClass& c : sn.classes)
        if (std::isinf(c.population)) return false;
    if (sn.nstations != 2) return false;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const SchedStrategy s = sn.stations[i].sched;
        if (s != SchedStrategy::PAS && s != SchedStrategy::OI) return false;
        if (!sn.stations[i].svc_rate_fun) return false;
    }
    return true;
}

/**
 * Port of `solver_nc_oi_analyzer.m`.
 *
 * @param sn  the refreshed struct
 * @param opt solver controls; unused, the analyzer is exact and has no tuning
 */
template <class T>
NcSolution<T> solver_nc_oi_analyzer(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    using qn::SchedStrategy;
    (void)opt;
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        throw UnsupportedError(
            "solver_nc_oi_analyzer: the BCMP weight table is formed as exp(lgamma(...)) and lG as "
            "log(G); this backend has no transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, K = sn.nclasses;

        // The OI rank rates are indexed by RAW class, so a chain that merges
        // several classes has no rate to evaluate.
        for (std::size_t c = 0; c < sn.nchains; ++c)
            if (sn.inchain[c].size() > 1)
                throw UnsupportedError(
                    "solver_nc_oi: requires one class per chain (no class switching)");
        std::vector<int> N(K, 0);
        for (std::size_t r = 0; r < K; ++r) {
            if (std::isinf(sn.classes[r].population))
                throw UnsupportedError("solver_nc_oi: requires a closed queueing network");
            N[r] = static_cast<int>(std::llround(sn.classes[r].population));
        }

        std::vector<bool> isOI(M, false), isINF(M, false), isQ(M, false);
        for (std::size_t i = 0; i < M; ++i) {
            const SchedStrategy s = sn.stations[i].sched;
            if (s == SchedStrategy::INF) {
                isINF[i] = true;
            } else if (s == SchedStrategy::PAS || s == SchedStrategy::OI) {
                if (!qn::station_swap_graph_is_zero(sn, i + 1))
                    throw UnsupportedError(
                        "solver_nc_oi: supports OI stations only (PAS with a non-empty swap graph "
                        "is not order-independent)");
                isOI[i] = true;
                if (!sn.stations[i].svc_rate_fun)
                    throw UnsupportedError(
                        "solver_nc_oi: an OI station has no service rate function; set it via "
                        "setServiceRateFunction");
            } else if (s == SchedStrategy::PS || s == SchedStrategy::LCFSPR ||
                       s == SchedStrategy::FCFS || s == SchedStrategy::SIRO) {
                isQ[i] = true;
                if ((s == SchedStrategy::FCFS || s == SchedStrategy::SIRO) &&
                    detail::class_dependent_fcfs_rate(sn, i))
                    throw UnsupportedError(
                        "solver_nc_oi: a station has class-dependent FCFS/SIRO rates and is not "
                        "product form; class-independent rates are required");
            } else {
                throw UnsupportedError(
                    "solver_nc_oi: supports only INF (delay), OI, PS, LCFS-PR, SIRO and "
                    "class-independent FCFS stations");
            }
        }

        const Matrix<T> V = detail::oi_visits(sn);

        // Delays aggregate into Z; every other BCMP queue keeps its own demand.
        std::vector<T> Z(K, zero);
        Matrix<T> D(M, K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                const double mu = num_traits<T>::to_double(sn.rates(i, r));
                if (!std::isfinite(mu) || mu == 0.0) continue;
                const T st = T(V(i, r) / sn.rates(i, r));
                if (isINF[i]) Z[r] += st;
                if (isQ[i]) D(i, r) = st;
            }

        for (std::size_t i = 0; i < M; ++i)
            if (isOI[i])
                for (std::size_t r = 0; r < K; ++r)
                    if (N[r] > 0 &&
                        std::fabs(num_traits<T>::to_double(V(i, r)) - 1.0) > 1e-9)
                        throw UnsupportedError(
                            "solver_nc_oi: requires unit per-class visits at every OI station");

        std::vector<std::size_t> oiList, qList;
        for (std::size_t i = 0; i < M; ++i) {
            if (isOI[i]) oiList.push_back(i + 1);
            if (isQ[i]) qList.push_back(i + 1);
        }
        std::vector<pfqn::OiRate<T>> rates;
        for (std::size_t m = 0; m < oiList.size(); ++m) {
            const std::function<T(const std::vector<std::size_t>&)> f =
                sn.stations[oiList[m] - 1].svc_rate_fun;
            rates.push_back(
                [f](const std::vector<int>& n) { return f(detail::oi_microstate(n)); });
        }

        std::vector<std::size_t> shp, stride;
        std::size_t total = 0;
        detail::oi_lattice(N, shp, stride, total);

        // The core table over the lattice: OI stations plus the aggregated delay.
        std::vector<T> Gfull(total, zero);
        for (std::size_t i = 0; i < total; ++i)
            Gfull[i] = pfqn::pfqn_ncoi(Z, detail::oi_sub(i, shp), rates).G;

        for (std::size_t i : qList) {
            std::vector<T> Dq(K, zero);
            for (std::size_t r = 0; r < K; ++r) Dq[r] = D(i - 1, r);
            const std::vector<T> Wq =
                detail::oi_ld_table(Dq, sn.stations[i - 1].nservers, shp, total);
            Gfull = detail::oi_conv(Gfull, Wq, shp, stride, total);
        }

        const T G = Gfull[total - 1];
        out.sol.lG = num_traits<T>::log_as_double(G);

        std::vector<T> X(K, zero);
        for (std::size_t r = 0; r < K; ++r)
            if (N[r] > 0) {
                std::vector<int> Nr = N;
                --Nr[r];
                X[r] = T(Gfull[detail::oi_idx(Nr, stride)] / G);
            }

        Matrix<T> Q(M, K, zero), Tp(M, K, zero), R(M, K, zero), U(M, K, zero);
        const T one = num_traits<T>::from_int(1);
        for (std::size_t m = 0; m < oiList.size(); ++m) {
            const std::size_t i = oiList[m];
            const std::vector<T> Phi = detail::oi_phi<T>(rates[m], N);
            for (std::size_t r = 0; r < K; ++r) {
                if (N[r] == 0) continue;
                const pfqn::OiFncResult<T> fr = pfqn::pfqn_oi_fnc<T>(
                    Phi, N, [r](const std::vector<int>& n) {
                        return num_traits<T>::from_int(n[r]);
                    });
                Q(i - 1, r) =
                    T(detail::oi_fnc_mean(fr.Psi, Gfull, shp, stride, total) / G - one);
            }
        }
        for (std::size_t i : qList) {
            std::vector<T> Dq(K, zero);
            for (std::size_t r = 0; r < K; ++r) Dq[r] = D(i - 1, r);
            const std::vector<T> Wq =
                detail::oi_ld_table(Dq, sn.stations[i - 1].nservers, shp, total);
            for (std::size_t r = 0; r < K; ++r) {
                if (N[r] == 0) continue;
                const pfqn::OiFncResult<T> fr = pfqn::pfqn_oi_fnc<T>(
                    Wq, N, [r](const std::vector<int>& n) {
                        return num_traits<T>::from_int(n[r]);
                    });
                Q(i - 1, r) =
                    T(detail::oi_fnc_mean(fr.Psi, Gfull, shp, stride, total) / G - one);
            }
        }
        for (std::size_t i = 0; i < M; ++i)
            if (isINF[i])
                for (std::size_t r = 0; r < K; ++r) {
                    const double mu = num_traits<T>::to_double(sn.rates(i, r));
                    if (!std::isfinite(mu) || mu == 0.0) continue;
                    Q(i, r) = T(X[r] * V(i, r) / sn.rates(i, r));
                }

        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) Tp(i, r) = T(X[r] * V(i, r));
        for (std::size_t i = 0; i < M; ++i)
            if (isINF[i])
                for (std::size_t r = 0; r < K; ++r) U(i, r) = Q(i, r);  // INF convention
        for (std::size_t i : qList) {
            const double c_raw = sn.stations[i - 1].nservers;
            const double c = (std::isfinite(c_raw) && c_raw > 0.0) ? c_raw : 1.0;
            for (std::size_t r = 0; r < K; ++r)
                U(i - 1, r) = T(X[r] * D(i - 1, r) / num_traits<T>::from_double(c));
        }
        for (std::size_t m = 0; m < oiList.size(); ++m) {
            // IN-SERVICE utilization E[sir_r]/c, the exact CTMC / LDES
            // convention: sir_r counts the class-r jobs receiving a strictly
            // positive rank rate, which is a function of the count vector, so
            // its mean is read off the same functional-server identity. It
            // coincides with the offered-load form T/mu(e_r)/c only when a job
            // engages a single server.
            const std::size_t i = oiList[m];
            const double s_raw = sn.stations[i - 1].nservers;
            const double s = (std::isfinite(s_raw) && s_raw > 0.0) ? s_raw : 1.0;
            const std::vector<T> Phi = detail::oi_phi<T>(rates[m], N);
            const pfqn::OiInsvcResult<T> ins = pfqn::pfqn_oi_insvc<T>(rates[m], N);
            for (std::size_t r = 0; r < K; ++r) {
                if (N[r] == 0) continue;
                const pfqn::OiFncResult<T> fr = pfqn::pfqn_oi_fnc<T>(
                    Phi, N, [&ins, &stride, r](const std::vector<int>& n) {
                        return ins.g(detail::oi_idx(n, stride), r);
                    });
                U(i - 1, r) =
                    T((detail::oi_fnc_mean(fr.Psi, Gfull, shp, stride, total) / G - one) /
                      num_traits<T>::from_double(s));
            }
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                if (Tp(i, r) > zero) R(i, r) = T(Q(i, r) / Tp(i, r));

        std::vector<T> C(K, zero);
        for (std::size_t r = 0; r < K; ++r)
            if (X[r] > zero) C[r] = T(num_traits<T>::from_int(N[r]) / X[r]);

        out.sol.Q = Q;
        out.sol.U = U;
        out.sol.R = R;
        out.sol.Tp = Tp;
        out.sol.X = X;
        out.sol.C = C;
        out.sol.iter = 1;
        out.sol.method = "oi";
        out.actualmethod = "oi";
        return out;
    }
}

/**
 * Port of `solver_nc_pas_is_analyzer.m`.
 *
 * @param sn  the refreshed struct
 * @param opt solver controls; `samples` and `seed` reach `pfqn_pas_is`
 */
template <class T>
NcSolution<T> solver_nc_pas_is_analyzer(const qn::NetworkStruct<T>& sn,
                                        const NcSolverOptions& opt) {
    using qn::SchedStrategy;
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc_pas_is_analyzer: the auto-normalized importance sampler draws from a "
            "continuous proposal and reports lG = log(G); this backend has no transcendental "
            "arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, K = sn.nclasses;

        for (std::size_t c = 0; c < sn.nchains; ++c)
            if (sn.inchain[c].size() > 1)
                throw UnsupportedError(
                    "solver_nc_pas_is: requires one class per chain (no class switching)");
        for (const qn::JobClass& c : sn.classes)
            if (std::isinf(c.population))
                throw UnsupportedError("solver_nc_pas_is: requires a closed queueing network");
        if (M != 2)
            throw UnsupportedError(
                "solver_nc_pas_is: models a two-station pass-and-swap tandem");
        std::vector<int> N(K, 0);
        for (std::size_t r = 0; r < K; ++r)
            N[r] = static_cast<int>(std::llround(sn.classes[r].population));

        std::vector<std::function<T(const std::vector<std::size_t>&)>> svc(M);
        std::vector<Matrix<T>> swapG(M);
        for (std::size_t i = 0; i < M; ++i) {
            const SchedStrategy s = sn.stations[i].sched;
            if (s != SchedStrategy::PAS && s != SchedStrategy::OI)
                throw UnsupportedError(
                    "solver_nc_pas_is: requires both stations to be OI/PAS");
            svc[i] = sn.stations[i].svc_rate_fun;
            if (!svc[i])
                throw UnsupportedError(
                    "solver_nc_pas_is: an OI/PAS station has no service rate function; set it via "
                    "setServiceRateFunction");
            swapG[i] = qn::station_swap_graph(sn, i + 1);
        }

        std::vector<pfqn::OiRateFun<T>> mu(M);
        for (std::size_t i = 0; i < M; ++i) mu[i] = detail::oi_rank_rate<T>(svc[i]);

        // The stored graph is the raw class-compatibility graph; the estimator
        // needs the global placement-order DAG that defines the recurrent
        // communicating class, derived from the P&S dynamics on the one-job-per-
        // class instance. Station 2 is the reversed suffix inside pfqn_pas_is,
        // so a single global order suffices.
        //
        // pas_swap2order hands its rate functions an ORDERED microstate (a list
        // of 1-based class indices), where pfqn_pas_is hands its own a support
        // indicator; both reach the same svcRateFun, so only the adapter differs.
        std::vector<pfqn::PasRateFun<T>> murates(M);
        for (std::size_t i = 0; i < M; ++i) {
            const std::function<T(const std::vector<std::size_t>&)> f = svc[i];
            murates[i] = [f](const std::vector<int>& c) {
                std::vector<std::size_t> mc(c.size());
                for (std::size_t a = 0; a < c.size(); ++a)
                    mc[a] = static_cast<std::size_t>(c[a]);
                return f(mc);
            };
        }
        const Matrix<T> H = pfqn::pas_swap2order<T>(swapG, murates, std::vector<int>(K, 1));
        Matrix<int> Hi(H.rows(), H.cols(), 0);
        for (std::size_t a = 0; a < H.rows(); ++a)
            for (std::size_t b = 0; b < H.cols(); ++b)
                Hi(a, b) = static_cast<int>(std::llround(num_traits<T>::to_double(H(a, b))));

        const Matrix<T> V = detail::oi_visits(sn);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                if (N[r] > 0 && std::fabs(num_traits<T>::to_double(V(i, r)) - 1.0) > 1e-9)
                    throw UnsupportedError(
                        "solver_nc_pas_is: requires unit per-class visits at both stations");

        const std::size_t samples = opt.samples > 0 ? opt.samples : 10000;
        const unsigned seed = opt.seed > 0 ? static_cast<unsigned>(opt.seed) : 23456u;

        pfqn::McRng g0(seed);
        const pfqn::PasIsResult<T> res = pfqn::pfqn_pas_is<T>(N, mu, Hi, samples, g0);
        out.sol.lG = num_traits<T>::log_as_double(res.G);

        Matrix<T> Q(M, K, zero), Tp(M, K, zero), R(M, K, zero), U(M, K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) Q(i, r) = res.Q(i, r);

        // Common random numbers across the N and N-e_r runs: each call restarts
        // the stream from the same seed, which is what makes the ratio far less
        // noisy than two independent estimates would be.
        std::vector<T> X(K, zero);
        for (std::size_t r = 0; r < K; ++r)
            if (N[r] > 0) {
                std::vector<int> Nr = N;
                --Nr[r];
                pfqn::McRng gr_rng(seed);
                // Only the constant is read here, so this run skips the
                // prefix-count coefficients: same stream, same G, none of the
                // O(ell R) per-sample bookkeeping behind the queue lengths.
                const pfqn::PasIsResult<T> gr =
                    pfqn::pfqn_pas_is<T>(Nr, mu, Hi, samples, gr_rng, false);
                if (res.G > zero) X[r] = T(gr.G / res.G);
            }

        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) Tp(i, r) = T(X[r] * V(i, r));
        for (std::size_t i = 0; i < M; ++i) {
            const double s_raw = sn.stations[i].nservers;
            const double s = (std::isfinite(s_raw) && s_raw > 0.0) ? s_raw : 1.0;
            for (std::size_t r = 0; r < K; ++r) {
                if (N[r] == 0) continue;
                std::vector<int> er(K, 0);
                er[r] = 1;
                const T muR = mu[i](er);  // rank rate with only class r present
                if (muR > zero)
                    U(i, r) = T(Tp(i, r) / muR / num_traits<T>::from_double(s));
            }
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                if (Tp(i, r) > zero) R(i, r) = T(Q(i, r) / Tp(i, r));

        std::vector<T> C(K, zero);
        for (std::size_t r = 0; r < K; ++r)
            if (X[r] > zero) C[r] = T(num_traits<T>::from_int(N[r]) / X[r]);

        out.sol.Q = Q;
        out.sol.U = U;
        out.sol.R = R;
        out.sol.Tp = Tp;
        out.sol.X = X;
        out.sol.C = C;
        out.sol.iter = 1;
        out.sol.method = "is";
        out.actualmethod = "is";
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_OI_H
