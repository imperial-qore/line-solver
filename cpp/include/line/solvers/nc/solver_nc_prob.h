/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_PROB_H
#define LINE_SOLVERS_NC_SOLVER_NC_PROB_H

/**
 * The state-probability half of the SolverNC class surface: ports of
 * `solver_nc_marg.m`, `solver_nc_margaggr.m`, `solver_nc_joint.m`,
 * `solver_nc_jointaggr.m` and `solver_nc_jointaggr_ld.m`, with the five
 * `@@SolverNC/getProb*` entry points on top of them.
 *
 * WHAT THIS ADDS THAT NOTHING ELSE IN THE TREE HAS. These are EXACT
 * product-form state probabilities. `solver_mva_prob.h` also answers
 * `getProbAggr`, but by its own account it fits a binomial to the means
 * (Schmidt 1997); here the probability is a ratio of normalizing constants and
 * is the model's own, to the last bit.
 *
 * THE IDENTITY THEY ALL USE. For a station i holding the per-class vector n_i,
 *
 *     Pr[n_i] = F_i(n_i) G_{-i}(N - n_i) / G(N)
 *
 * where F_i is the station's own balance function evaluated at n_i, G_{-i} is
 * the constant of the network with station i deleted, and G is the constant of
 * the whole model. Each factor is another `pfqn_ncld` call on a load-dependent
 * lattice, so the whole family is three or four constants per station and
 * nothing else.
 *
 * NO STATE PACKAGE: THE INPUT IS THE MARGINAL VECTOR. The reference reaches the
 * per-class counts through `State.toMarginal(sn, ist, state{isf})`, and
 * `getProbAggr` gets there by encoding the user's per-class vector with
 * `State.fromMarginal` first -- a round trip whose only product is the vector
 * the user already supplied. This port takes that vector directly. The
 * consequence is precise and is enforced rather than hidden: two branches of
 * `solver_nc_marg` read state a marginal does not carry, and both are refused
 * by name (see `solver_nc_prob`).
 *
 * ARITHMETIC. Every probability is a difference of logarithms of normalizing
 * constants, exponentiated once; a non-transcendental backend is refused by
 * name, as in `solver_nc.h`.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_jointmarg.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/api/pfqn/pfqn_procomom.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/solver_nc.h"
#include "line/solvers/nc/solver_ncld.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/**
 * A state, as this port expresses it: `nir[i][r]` jobs of class r at station i.
 *
 * This is `State.toMarginal`'s second output and `State.fromMarginal`'s input,
 * i.e. the only part of the reference's state encoding these analyzers use. A
 * NEGATIVE entry is the reference's "ignore this station" flag and is honoured.
 */
using MarginalState = std::vector<std::vector<int>>;

namespace detail {

/** The load-dependent lattice `mu` the probability analyzers all build. */
template <class T>
Matrix<T> prob_mu(const qn::NetworkStruct<T>& sn, std::size_t Ntot) {
    const std::size_t M = sn.nstations;
    const std::size_t w = std::max<std::size_t>(1, Ntot);
    Matrix<T> mu(M, w, num_traits<T>::from_int(1));
    for (std::size_t i = 0; i < M; ++i) {
        const double S = sn.stations[i].nservers;
        for (std::size_t n = 1; n <= w; ++n)
            mu(i, n - 1) = num_traits<T>::from_double(
                std::isinf(S) ? static_cast<double>(n)
                              : std::min<double>(static_cast<double>(n), S));
    }
    return mu;
}

/** `nivec * sn.chains'`: the per-chain totals of a per-class vector. */
template <class T>
std::vector<int> to_chain(const qn::NetworkStruct<T>& sn, const std::vector<int>& nir) {
    std::vector<int> nc(sn.nchains, 0);
    for (std::size_t c = 0; c < sn.nchains; ++c)
        for (std::size_t k : sn.inchain[c]) nc[c] += nir[k - 1];
    return nc;
}

/** One row of a matrix, as the 1 x R matrix `pfqn_ncld` wants. */
template <class T>
Matrix<T> row_of(const Matrix<T>& A, std::size_t i) {
    Matrix<T> r(1, A.cols(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < A.cols(); ++j) r(0, j) = A(i, j);
    return r;
}

/** Every row but one (MATLAB's `A(setdiff(1:M,i),:)`). */
template <class T>
Matrix<T> drop_row(const Matrix<T>& A, std::size_t i) {
    Matrix<T> r(A.rows() - 1, A.cols(), num_traits<T>::from_int(0));
    std::size_t o = 0;
    for (std::size_t k = 0; k < A.rows(); ++k) {
        if (k == i) continue;
        for (std::size_t j = 0; j < A.cols(); ++j) r(o, j) = A(k, j);
        ++o;
    }
    return r;
}

/** `sn.rates` reciprocated, zero where the pair is disabled. */
template <class T>
Matrix<T> service_times(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> ST(sn.nstations, sn.nclasses, zero);
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            if (!sn.disabled[i][r] && sn.rates(i, r) != zero) ST(i, r) = T(one / sn.rates(i, r));
    return ST;
}

/**
 * `V(i,k)` as the probability analyzers build it: the visit of class k at
 * station i within ITS OWN chain, unnormalized.
 *
 * The reference indexes `sn.visits{c}(ist,k)` with a STATION index while
 * `sn.visits` is indexed by stateful node; the two coincide on every model NC
 * accepts today, since a Cache is refused and no other node is stateful. This
 * port uses the stateful index, which is what the field is.
 */
template <class T>
Matrix<T> class_visits(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> V(sn.nstations, sn.nclasses, zero);
    for (std::size_t c = 0; c < sn.nchains; ++c)
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            const std::size_t sf = sn.stateful_of_station(i + 1) - 1;
            for (std::size_t k : sn.inchain[c]) V(i, k - 1) = sn.visits[c](sf, k - 1);
        }
    return V;
}

/** The total closed population, refusing an open model the way these do. */
template <class T>
std::size_t closed_total(const qn::NetworkStruct<T>& sn, const char* who) {
    double t = 0.0;
    for (const qn::JobClass& c : sn.classes) {
        if (std::isinf(c.population))
            throw UnsupportedError(std::string(who) +
                                   ": the state probability is defined on a CLOSED network, and "
                                   "this model has an open class");
        t += c.population;
    }
    return static_cast<std::size_t>(std::llround(t));
}

/** Validate a marginal against the model, so a bad index is not a wrong number. */
template <class T>
void check_marginal(const qn::NetworkStruct<T>& sn, const MarginalState& nir, const char* who) {
    if (nir.size() != sn.nstations)
        throw InputError(std::string(who) + ": the marginal state has " +
                         std::to_string(nir.size()) + " stations, the model has " +
                         std::to_string(sn.nstations));
    for (const std::vector<int>& row : nir)
        if (row.size() != sn.nclasses)
            throw InputError(std::string(who) +
                             ": every station's marginal must give one count per class");
}

}  // namespace detail

/** What the marginal analyzers return: one probability per station. */
template <class T>
struct NcMargResult {
    std::vector<T> P;     ///< (M) probability that station i holds its given vector
    std::vector<T> logP;  ///< (M) the same, in logs
    double lG = 0.0;      ///< the log normalizing constant that normalized them
};

/**
 * Port of `solver_nc_margaggr.m`.
 *
 * The purely AGGREGATE marginal: the station's balance function evaluated at
 * the per-class vector, times the constant of the network without it. It reads
 * nothing but the marginal, so it carries no discipline restriction at all --
 * which is why `getProbAggr`, `getProbMarg` and `getProbSysAggr` are
 * unrestricted while `getProb` is not.
 *
 * @param sn  the refreshed struct
 * @param opt solver controls
 * @param nir the state; a station whose row has a NEGATIVE entry is skipped and
 *            reported as probability zero, which is the reference's flag
 * @param lG  a precomputed log normalizing constant; NaN to compute one
 */
template <class T>
NcMargResult<T> solver_nc_margaggr(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                                   const MarginalState& nir, double lG) {
    NcMargResult<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn; (void)opt; (void)nir; (void)lG;
        throw UnsupportedError(
            "solver_nc_margaggr: the state probability is exp(lF_i + lG_{-i} - lG), a difference "
            "of logarithms of normalizing constants, and needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, K = sn.nclasses;
        detail::check_marginal(sn, nir, "solver_nc_margaggr");
        const std::size_t Ntot = detail::closed_total(sn, "solver_nc_margaggr");

        const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        const Matrix<T> ST = detail::service_times(sn);
        const Matrix<T> V = detail::class_visits(sn);
        const Matrix<T> mu = detail::prob_mu(sn, Ntot);
        std::vector<int> Nchain(sn.nchains, 0);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            Nchain[c] = static_cast<int>(std::llround(d.Nchain[c]));

        const Matrix<T> Zc(1, sn.nchains, zero);
        const Matrix<T> Zk(1, K, zero);
        const pfqn::NcldMethod pm = detail::ncld_pfqn_method(opt.method);
        pfqn::NcOptions nopt;
        nopt.samples = opt.samples;
        nopt.seed = opt.seed;
        nopt.tol = opt.tol;
        const T atol = num_traits<T>::from_double(opt.tol);
        if (std::isnan(lG))
            lG = pfqn::pfqn_ncld(d.Lchain, Nchain, Zc, mu, pm, atol, nopt).lG;
        out.lG = lG;

        out.P.assign(M, zero);
        out.logP.assign(M, zero);
        for (std::size_t i = 0; i < M; ++i) {
            bool ignore = false;
            for (int v : nir[i])
                if (v < 0) ignore = true;
            if (ignore) continue;  // MATLAB sets NaN here and then Pr(isnan)=0
            const std::vector<int> nc = detail::to_chain(sn, nir[i]);
            std::vector<int> Nrest(sn.nchains, 0);
            for (std::size_t c = 0; c < sn.nchains; ++c) Nrest[c] = Nchain[c] - nc[c];
            const double lG_minus_i =
                M > 1 ? pfqn::pfqn_ncld(detail::drop_row(d.Lchain, i), Nrest, Zc,
                                        detail::drop_row(mu, i), pm, atol, nopt)
                            .lG
                      : 0.0;
            Matrix<T> Fi(1, K, zero);
            for (std::size_t r = 0; r < K; ++r) Fi(0, r) = T(ST(i, r) * V(i, r));
            const double lF_i =
                pfqn::pfqn_ncld(Fi, nir[i], Zk, detail::row_of(mu, i), pm, atol, nopt).lG;
            const double lp = lF_i + lG_minus_i - lG;
            out.logP[i] = num_traits<T>::from_double(lp);
            out.P[i] = num_traits<T>::from_double(std::exp(lp));
        }
        return out;
    }
}

/**
 * Port of `solver_nc_marg.m`: the DETAILED marginal, which weighs the station's
 * internal arrangement and therefore depends on its discipline.
 *
 * TWO BRANCHES ARE REFUSED BY NAME because a per-class marginal cannot carry
 * what they read, and answering them from the default arrangement would be a
 * fabricated number:
 *
 *   SIRO wants the CLASS OF THE JOB IN SERVICE (`sivec`), whose term is
 *        log(n_ci / sum n). A marginal says how many jobs of each class are
 *        present, not which one holds the server.
 *   PS and INF want the PHASE-LEVEL occupancy (`kirvec`) when service is not
 *        exponential. Under exponential service kirvec IS the marginal, and
 *        that case is computed exactly.
 *
 * FCFS additionally carries the reference's own preconditions -- exponential
 * service, and identical mean service time across the classes -- without which
 * the station is not product-form and the reference errors out.
 */
template <class T>
NcMargResult<T> solver_nc_marg(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                               const MarginalState& nir, double lG) {
    NcMargResult<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn; (void)opt; (void)nir; (void)lG;
        throw UnsupportedError(
            "solver_nc_marg: the state probability is exp(lF_i + lG_{-i} - lG), a difference of "
            "logarithms of normalizing constants, and needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        const std::size_t M = sn.nstations, K = sn.nclasses;
        detail::check_marginal(sn, nir, "solver_nc_marg");
        const std::size_t Ntot = detail::closed_total(sn, "solver_nc_marg");

        const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        const Matrix<T> ST = detail::service_times(sn);
        const Matrix<T> V = detail::class_visits(sn);
        const Matrix<T> mu = detail::prob_mu(sn, Ntot);
        std::vector<int> Nchain(sn.nchains, 0);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            Nchain[c] = static_cast<int>(std::llround(d.Nchain[c]));

        const Matrix<T> Zc(1, sn.nchains, zero);
        const Matrix<T> Zk(1, K, zero);
        const pfqn::NcldMethod pm = detail::ncld_pfqn_method(opt.method);
        pfqn::NcOptions nopt;
        nopt.samples = opt.samples;
        nopt.seed = opt.seed;
        nopt.tol = opt.tol;
        const T atol = num_traits<T>::from_double(opt.tol);
        if (std::isnan(lG)) lG = pfqn::pfqn_ncld(d.Lchain, Nchain, Zc, mu, pm, atol, nopt).lG;
        out.lG = lG;

        // A station is exponential in class r when its service law has one phase.
        const auto is_exponential = [&](std::size_t i, std::size_t r) {
            if (sn.disabled[i][r]) return true;
            return sn.service[i][r].D0.rows() <= 1;
        };

        out.P.assign(M, zero);
        out.logP.assign(M, zero);
        for (std::size_t i = 0; i < M; ++i) {
            bool ignore = false;
            for (int v : nir[i])
                if (v < 0) ignore = true;
            if (ignore) continue;
            const std::vector<int> nc = detail::to_chain(sn, nir[i]);
            std::vector<int> Nrest(sn.nchains, 0);
            for (std::size_t c = 0; c < sn.nchains; ++c) Nrest[c] = Nchain[c] - nc[c];
            const double lG_minus_i =
                M > 1 ? pfqn::pfqn_ncld(detail::drop_row(d.Lchain, i), Nrest, Zc,
                                        detail::drop_row(mu, i), pm, atol, nopt)
                            .lG
                      : 0.0;

            long ntot_i = 0;
            for (int v : nir[i]) ntot_i += v;
            // sum_{n=1}^{|n_i|} log mu_i(n), the load-dependent denominator
            double lmu = 0.0;
            for (long n = 1; n <= ntot_i && n <= static_cast<long>(mu.cols()); ++n)
                lmu += std::log(num_traits<T>::to_double(mu(i, static_cast<std::size_t>(n - 1))));

            double lF_i = 0.0;
            const qn::SchedStrategy sc = sn.stations[i].sched;
            if (sc == qn::SchedStrategy::FCFS) {
                double stmax = 0.0;
                for (std::size_t r = 0; r < K; ++r) {
                    if (sn.disabled[i][r]) continue;
                    if (!is_exponential(i, r))
                        throw UnsupportedError(
                            "solver_nc_marg: the product-form state probability requires "
                            "exponential service times at FCFS nodes, and this station's class " +
                            std::to_string(r + 1) + " is not exponential");
                    stmax = std::max(stmax, num_traits<T>::to_double(ST(i, r)));
                }
                for (std::size_t r = 0; r < K; ++r) {
                    if (sn.disabled[i][r] || nir[i][r] == 0) continue;
                    if (std::fabs(num_traits<T>::to_double(ST(i, r)) - stmax) >
                        GlobalConstants::FineTol)
                        throw UnsupportedError(
                            "solver_nc_marg: the product-form state probability requires "
                            "identical service times across classes at FCFS nodes, and this "
                            "station's class " + std::to_string(r + 1) + " differs");
                }
                if (ntot_i > 0) {
                    // REFERENCE DEFECT, corrected here: the reference writes
                    // `sum(nirvec .* log(V(ist,r)))` with `r` left over from the
                    // validation loop above it, so EVERY class is weighted by the
                    // LAST class's visit ratio. The intended term is each class's
                    // own visit, which is what the identity Pr[n] ~ prod_r V_ir^{n_ir}
                    // requires and what every other branch of this file uses.
                    for (std::size_t r = 0; r < K; ++r) {
                        if (nir[i][r] == 0) continue;
                        const double v = num_traits<T>::to_double(V(i, r));
                        if (!(v > 0.0))
                            throw NumericError(
                                "solver_nc_marg: class " + std::to_string(r + 1) +
                                " holds jobs at a station it never visits");
                        lF_i += static_cast<double>(nir[i][r]) * std::log(v);
                    }
                    lF_i -= lmu;
                }
            } else if (sc == qn::SchedStrategy::SIRO) {
                throw UnsupportedError(
                    "solver_nc_marg: the SIRO branch weighs the state by log(n_ci / sum n), which "
                    "needs the CLASS OF THE JOB IN SERVICE; a per-class marginal does not carry "
                    "it. Use getProbAggr, whose aggregate marginal has no such dependency");
            } else if (sc == qn::SchedStrategy::PS || sc == qn::SchedStrategy::INF) {
                for (std::size_t r = 0; r < K; ++r) {
                    if (sn.disabled[i][r]) continue;
                    if (!is_exponential(i, r))
                        throw UnsupportedError(
                            "solver_nc_marg: a non-exponential service law at a " +
                            std::string(sc == qn::SchedStrategy::PS ? "PS" : "delay") +
                            " station makes the balance function depend on the PHASE-LEVEL "
                            "occupancy, which a per-class marginal does not carry");
                    if (nir[i][r] == 0) continue;
                    const double w = num_traits<T>::to_double(T(V(i, r) * ST(i, r)));
                    if (!(w > 0.0))
                        throw NumericError(
                            "solver_nc_marg: class " + std::to_string(r + 1) +
                            " holds jobs at a station whose demand for it is zero");
                    lF_i += static_cast<double>(nir[i][r]) * std::log(w);
                    for (int q = 2; q <= nir[i][r]; ++q) lF_i -= std::log(static_cast<double>(q));
                }
                for (long q = 2; q <= ntot_i; ++q) lF_i += std::log(static_cast<double>(q));
                lF_i -= lmu;
            }
            // Any other discipline leaves lF_i at zero, as the reference's
            // switch does when no case matches.
            (void)one;

            const double lp = lF_i + lG_minus_i - lG;
            out.logP[i] = num_traits<T>::from_double(lp);
            out.P[i] = num_traits<T>::from_double(std::exp(lp));
        }
        return out;
    }
}

/**
 * Port of `solver_nc_joint.m`: the probability of the WHOLE system state.
 *
 * The per-station factor is the chain-level balance function corrected by the
 * class-within-chain split `lg0_i - lG0_i`, which is what turns a chain-level
 * constant into a class-level one.
 */
template <class T>
T solver_nc_joint(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                  const MarginalState& nir, double* lG_out) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn; (void)opt; (void)nir; (void)lG_out;
        throw UnsupportedError(
            "solver_nc_joint: the joint state probability is a difference of logarithms of "
            "normalizing constants and needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, K = sn.nclasses;
        detail::check_marginal(sn, nir, "solver_nc_joint");
        const std::size_t Ntot = detail::closed_total(sn, "solver_nc_joint");

        const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        const Matrix<T> ST = detail::service_times(sn);
        const Matrix<T> mu = detail::prob_mu(sn, Ntot);
        std::vector<int> Nchain(sn.nchains, 0);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            Nchain[c] = static_cast<int>(std::llround(d.Nchain[c]));

        const Matrix<T> Zc(1, sn.nchains, zero);
        const Matrix<T> Zk(1, K, zero);
        const pfqn::NcldMethod pm = detail::ncld_pfqn_method(opt.method);
        pfqn::NcOptions nopt;
        nopt.samples = opt.samples;
        nopt.seed = opt.seed;
        nopt.tol = opt.tol;
        const T atol = num_traits<T>::from_double(opt.tol);
        const double lG = pfqn::pfqn_ncld(d.Lchain, Nchain, Zc, mu, pm, atol, nopt).lG;
        if (lG_out != nullptr) *lG_out = lG;

        double lPr = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            const std::vector<int> nc = detail::to_chain(sn, nir[i]);
            const Matrix<T> mui = detail::row_of(mu, i);
            const double lF_i =
                pfqn::pfqn_ncld(detail::row_of(d.Lchain, i), nc, Zc, mui, pm, atol, nopt).lG;
            Matrix<T> g0(1, K, zero);
            for (std::size_t r = 0; r < K; ++r) g0(0, r) = T(ST(i, r) * d.alpha(i, r));
            const double lg0_i = pfqn::pfqn_ncld(g0, nir[i], Zk, mui, pm, atol, nopt).lG;
            const double lG0_i =
                pfqn::pfqn_ncld(detail::row_of(d.STchain, i), nc, Zc, mui, pm, atol, nopt).lG;
            lPr += lF_i + (lg0_i - lG0_i);
        }
        return num_traits<T>::from_double(std::exp(lPr - lG));
    }
}

/**
 * Port of `solver_nc_jointaggr.m`: the aggregate joint.
 *
 * The reference takes its constant from `pfqn_ncld` under `method='exact'` and
 * from `solver_nc` otherwise, noting in a comment that the second is "unclear
 * ... as it doesn't consider the transformation to ld model". Both paths are
 * reproduced, because the choice changes the answer and the caller's method is
 * what selects it.
 */
template <class T>
T solver_nc_jointaggr(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                      const MarginalState& nir, double* lG_out) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn; (void)opt; (void)nir; (void)lG_out;
        throw UnsupportedError(
            "solver_nc_jointaggr: the joint state probability is a difference of logarithms of "
            "normalizing constants and needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, K = sn.nclasses;
        detail::check_marginal(sn, nir, "solver_nc_jointaggr");
        const std::size_t Ntot = detail::closed_total(sn, "solver_nc_jointaggr");

        const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        const Matrix<T> mu = detail::prob_mu(sn, Ntot);
        std::vector<int> Nchain(sn.nchains, 0);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            Nchain[c] = static_cast<int>(std::llround(d.Nchain[c]));

        const Matrix<T> Zc(1, sn.nchains, zero);
        const Matrix<T> Zk(1, K, zero);
        const pfqn::NcldMethod pm = detail::ncld_pfqn_method(opt.method);
        pfqn::NcOptions nopt;
        nopt.samples = opt.samples;
        nopt.seed = opt.seed;
        nopt.tol = opt.tol;
        const T atol = num_traits<T>::from_double(opt.tol);

        double lG;
        Matrix<T> ST;
        if (opt.method == "exact") {
            lG = pfqn::pfqn_ncld(d.Lchain, Nchain, Zc, mu, pm, atol, nopt).lG;
            ST = detail::service_times(sn);
        } else {
            const NcSolution<T> s = solver_nc(sn, opt);
            lG = s.sol.lG;
            ST = s.STeff;
        }
        if (lG_out != nullptr) *lG_out = lG;

        // V is cellsum(sn.visits) here, summed over chains, not the per-chain
        // class visit the detailed marginal uses.
        Matrix<T> V(M, K, zero);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            for (std::size_t i = 0; i < M; ++i) {
                const std::size_t sf = sn.stateful_of_station(i + 1) - 1;
                for (std::size_t r = 0; r < K; ++r) V(i, r) = T(V(i, r) + sn.visits[c](sf, r));
            }

        double lPr = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            const std::vector<int> nc = detail::to_chain(sn, nir[i]);
            bool any = false;
            for (int v : nc)
                if (v > 0) any = true;
            if (!any) continue;
            Matrix<T> Fi(1, K, zero);
            for (std::size_t r = 0; r < K; ++r) Fi(0, r) = T(ST(i, r) * V(i, r));
            lPr += pfqn::pfqn_ncld(Fi, nir[i], Zk, detail::row_of(mu, i), pm, atol, nopt).lG;
        }
        return num_traits<T>::from_double(std::exp(lPr - lG));
    }
}

/**
 * Port of `solver_nc_jointaggr_ld.m`: the load-dependent joint.
 *
 * Identical to `solver_nc_joint` except that the chain demand comes straight
 * from `sn_get_demands_chain` rather than being rebuilt, which is what makes it
 * the load-dependent variant in the reference.
 */
template <class T>
T solver_nc_jointaggr_ld(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                         const MarginalState& nir, double* lG_out) {
    return solver_nc_joint(sn, opt, nir, lG_out);
}

// ---------------------------------------------------------------------------
// The @@SolverNC class surface
// ---------------------------------------------------------------------------

/**
 * Port of `@@SolverNC/getProb.m`: the DETAILED state probability at one station.
 *
 * RETURNS THE LOG PROBABILITY, NOT THE PROBABILITY, DESPITE THE NAME. The value
 * is therefore NEGATIVE and is not in [0,1]; `exp()` recovers the probability.
 * On Delay(1)+PS(2) at N=3 with two jobs at the queue this returns
 * -1.15267950993839, and the probability is exp of it, 0.31578947368.
 *
 * This is a DELIBERATE reproduction of the reference, not an oversight.
 * `solver_nc_marg.m` returns `lPr` as its first output and
 * `@@SolverNC/getProb.m` passes it through under the name `Pnir` without
 * exponentiating. The port originally returned the probability and exposed the
 * log separately; the user ruled for strict bug-for-bug parity with MATLAB over
 * the safer API, and this is that decision. See register row N8, which records
 * the argument on both sides.
 *
 * `getProbAggr` is unaffected and DOES return a probability, so the two
 * accessors disagree in kind -- another reason the reference behaviour is
 * surprising rather than merely unusual. Use `solver_nc_marg` directly for a
 * result carrying both `P` and `logP`.
 *
 * @param ist 1-based station index
 * @param nir the per-class occupancy at that station; other stations are left
 *            at the reference's "ignore" flag, since the detailed marginal is
 *            reported per station and only this one is asked for
 * @param sn the refreshed network struct
 * @param opt SolverNC's options
 * @return the LOG of the probability
 */
template <class T>
T solver_nc_getprob(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt, std::size_t ist,
                    const std::vector<int>& nir) {
    if (ist == 0 || ist > sn.nstations)
        throw InputError("getProb: station index out of range");
    MarginalState st(sn.nstations, std::vector<int>(sn.nclasses, -1));
    st[ist - 1] = nir;
    return solver_nc_marg(sn, opt, st,
                          std::numeric_limits<double>::quiet_NaN())
        .logP[ist - 1];
}

/** Port of `@@SolverNC/getProbAggr.m`: the AGGREGATE probability at one station. */
template <class T>
T solver_nc_getprob_aggr(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                         std::size_t ist, const std::vector<int>& nir, double lG) {
    if (ist == 0 || ist > sn.nstations)
        throw InputError("getProbAggr: station index out of range");
    MarginalState st(sn.nstations, std::vector<int>(sn.nclasses, -1));
    st[ist - 1] = nir;
    return solver_nc_margaggr(sn, opt, st, lG).P[ist - 1];
}

/** The marginal queue-length distribution and its logarithm. */
template <class T>
struct NcQueueLengthDist {
    std::vector<T> P;     ///< P[n] = Pr[n jobs at the station], n = 0..sum(N)
    std::vector<T> logP;
};

/**
 * Port of `@@SolverNC/getProbMarg.m`: the TOTAL queue-length distribution.
 *
 * Two routes, as in the reference. Under `method='comom'` the whole vector
 * comes from one `pfqn_procomom` solve; otherwise every total n is written as a
 * sum over the per-class partitions of n and each partition goes through the
 * aggregate marginal. The normalizing constant is computed ONCE and reused
 * across the partitions, which is the reference's caching and matters here:
 * the enumeration is otherwise quadratic in constants.
 */
template <class T>
NcQueueLengthDist<T> solver_nc_getprob_marg(const qn::NetworkStruct<T>& sn,
                                            const NcSolverOptions& opt, std::size_t ist) {
    NcQueueLengthDist<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn; (void)opt; (void)ist;
        throw UnsupportedError(
            "getProbMarg: the queue-length distribution is assembled from exponentiated "
            "log-constants and needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        if (ist == 0 || ist > sn.nstations)
            throw InputError("getProbMarg: station index out of range");
        const std::size_t K = sn.nclasses;
        const std::size_t Ntot = detail::closed_total(sn, "getProbMarg");

        std::vector<int> N(K, 0);
        for (std::size_t r = 0; r < K; ++r)
            N[r] = static_cast<int>(std::llround(sn.classes[r].population));

        out.P.assign(Ntot + 1, zero);
        out.logP.assign(Ntot + 1, num_traits<T>::from_double(
                                      -std::numeric_limits<double>::infinity()));

        if (opt.method == "comom") {
            // The fast path: pfqn_procomom returns the whole per-station
            // queue-length vector in one solve, on the Seidmann-reduced model.
            const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
            const std::size_t M = sn.nstations, C = sn.nchains;
            Matrix<T> Lms(M, C, zero);
            Matrix<T> Ztot(1, C, zero);
            std::vector<std::size_t> queueStations;
            for (std::size_t i = 0; i < M; ++i) {
                const double S = sn.stations[i].nservers;
                if (std::isinf(S)) {
                    for (std::size_t c = 0; c < C; ++c) Ztot(0, c) += d.Lchain(i, c);
                } else {
                    queueStations.push_back(i);
                    const T cs = num_traits<T>::from_double(S);
                    for (std::size_t c = 0; c < C; ++c) {
                        Lms(i, c) = T(d.Lchain(i, c) / cs);
                        Ztot(0, c) += T(d.Lchain(i, c) * num_traits<T>::from_double(S - 1.0) / cs);
                    }
                }
            }
            const auto pos = std::find(queueStations.begin(), queueStations.end(), ist - 1);
            if (pos == queueStations.end()) {
                // A DELAY station has no row in the procomom result. The
                // reference warns and FALLS BACK to the enumeration under
                // method='default' rather than refusing, so this does too.
                NcSolverOptions fallback = opt;
                fallback.method = "default";
                return solver_nc_getprob_marg(sn, fallback, ist);
            }
            Matrix<T> Lq(queueStations.size(), C, zero);
            for (std::size_t a = 0; a < queueStations.size(); ++a)
                for (std::size_t c = 0; c < C; ++c) Lq(a, c) = Lms(queueStations[a], c);
            std::vector<int> Nchain(C, 0);
            std::size_t sumNchain = 0;
            for (std::size_t c = 0; c < C; ++c) {
                Nchain[c] = static_cast<int>(std::llround(d.Nchain[c]));
                sumNchain += static_cast<std::size_t>(Nchain[c]);
            }
            std::vector<T> Zv(C, zero);
            for (std::size_t c = 0; c < C; ++c) Zv[c] = Ztot(0, c);
            const Matrix<T> Pr = pfqn::pfqn_procomom(Lq, Nchain, Zv).Pr;
            const std::size_t row = static_cast<std::size_t>(pos - queueStations.begin());
            const std::size_t len = std::min(sumNchain + 1, Ntot + 1);
            for (std::size_t n = 0; n < len && n < Pr.cols(); ++n) {
                out.P[n] = Pr(row, n);
                if (num_traits<T>::to_double(out.P[n]) > 0.0)
                    out.logP[n] = num_traits<T>::from_double(
                        std::log(num_traits<T>::to_double(out.P[n])));
            }
            return out;
        }

        // The enumeration. lG is computed once and handed to every call.
        const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        const Matrix<T> mu = detail::prob_mu(sn, Ntot);
        std::vector<int> Nchain(sn.nchains, 0);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            Nchain[c] = static_cast<int>(std::llround(d.Nchain[c]));
        const double lG =
            pfqn::pfqn_ncld(d.Lchain, Nchain, Matrix<T>(1, sn.nchains, zero), mu,
                            detail::ncld_pfqn_method(opt.method),
                            num_traits<T>::from_double(opt.tol))
                .lG;

        std::vector<int> part(K, 0);
        for (std::size_t n = 0; n <= Ntot; ++n) {
            double acc = 0.0;
            bool any = false;
            // enumerate the compositions of n into K parts with part_r <= N_r
            std::function<void(std::size_t, int)> rec = [&](std::size_t r, int left) {
                if (r + 1 == K) {
                    if (left > N[r]) return;
                    part[r] = left;
                    const T p = solver_nc_getprob_aggr(sn, opt, ist, part, lG);
                    const double v = num_traits<T>::to_double(p);
                    if (v > 0.0) {
                        acc += v;
                        any = true;
                    }
                    return;
                }
                const int hi = std::min(left, N[r]);
                for (int v = 0; v <= hi; ++v) {
                    part[r] = v;
                    rec(r + 1, left - v);
                }
            };
            if (K == 0) continue;
            rec(0, static_cast<int>(n));
            if (any) {
                out.P[n] = num_traits<T>::from_double(acc);
                out.logP[n] = num_traits<T>::from_double(std::log(acc));
            }
        }
        // THE REFERENCE RENORMALIZES, and so must this: `getProbMarg.m` divides
        // the whole vector by its sum whenever that sum misses one by more than
        // 1e-10. The gap is the normalizing constant's own error -- on the
        // default (`cub`) route lG is an approximation, and the per-partition
        // constants inherit it -- so without this the marginal law of a station
        // is not a distribution. The 1e-10 dead band is the reference's too: it
        // leaves a vector already summing to one bit-for-bit alone rather than
        // perturbing it by a division.
        double total = 0.0;
        for (std::size_t n = 0; n <= Ntot; ++n) total += num_traits<T>::to_double(out.P[n]);
        if (total > 0.0 && std::fabs(total - 1.0) > 1e-10) {
            const double ltotal = std::log(total);
            for (std::size_t n = 0; n <= Ntot; ++n) {
                out.P[n] = num_traits<T>::from_double(num_traits<T>::to_double(out.P[n]) / total);
                if (num_traits<T>::to_double(out.P[n]) > 0.0)
                    out.logP[n] =
                        num_traits<T>::from_double(num_traits<T>::to_double(out.logP[n]) - ltotal);
            }
        }
        return out;
    }
}

/** Port of `@@SolverNC/getProbSys.m`. */
template <class T>
T solver_nc_getprob_sys(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                        const MarginalState& nir) {
    return solver_nc_joint(sn, opt, nir, nullptr);
}

/** Port of `@@SolverNC/getProbSysAggr.m`. */
template <class T>
T solver_nc_getprob_sys_aggr(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                             const MarginalState& nir) {
    return solver_nc_jointaggr(sn, opt, nir, nullptr);
}

namespace detail {

/**
 * The permanent identity supplies one n_i! per queueing station and none per
 * infinite server. A multiserver or load-dependent station has neither, so it
 * is refused BY NAME rather than approximated.
 */
template <class T>
void check_jointmarg_supported(const qn::NetworkStruct<T>& sn) {
    for (const qn::JobClass& c : sn.classes)
        if (std::isinf(c.population))
            throw UnsupportedError(
                "solver_nc_jointmarg: getProbSysMarg requires a closed model: the joint law of "
                "the total queue lengths is not defined when a class has an infinite population");
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (!sn.stations[i].lldscaling.empty())
            throw UnsupportedError(
                "solver_nc_jointmarg: getProbSysMarg does not support load-dependent stations "
                "(sn.lldscaling is set at station " + std::to_string(i + 1) + "): the permanent "
                "identity supplies exactly one n_i! per queueing station");
        const double S = sn.stations[i].nservers;
        if (std::isfinite(S) && S > 1.0)
            throw UnsupportedError(
                "solver_nc_jointmarg: getProbSysMarg does not support the multiserver station " +
                std::to_string(i + 1) + " (" + std::to_string(static_cast<int>(S)) +
                " servers): the permanent identity supplies exactly one n_i! per queueing "
                "station");
    }
}

}  // namespace detail

/**
 * Joint probability that station i holds `nvec[i]` jobs IN TOTAL, all classes
 * summed out.
 *
 * This is NOT `solver_nc_jointaggr`, which fixes the per-class population of
 * every station: each state here is the SUM of jointaggr over the whole fibre
 * of per-class tables with these row sums, and that fibre grows
 * combinatorially. `pfqn_jointmarg` evaluates the sum in closed form as a
 * permanent of the demand matrix replicated once per job.
 *
 * @param nvec   (M) per-station total job counts
 * @param engine "exact" (default), "spm", "bethe", "heur", "huberlaw" or
 *               "adapart"; see pfqn_jointmarg for what each guarantees
 * @param lG_out receives the log normalizing constant when not null
 */
template <class T>
T solver_nc_jointmarg(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                      const std::vector<int>& nvec, const std::string& engine = "exact",
                      double* lG_out = nullptr) {
    const std::size_t M = sn.nstations;
    if (nvec.size() != M)
        throw InputError("solver_nc_jointmarg: the occupancy vector has " +
                         std::to_string(nvec.size()) + " entries but the model has " +
                         std::to_string(M) + " stations");
    detail::check_jointmarg_supported(sn);

    const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Lchain = d.Lchain;
    for (std::size_t i = 0; i < Lchain.rows(); ++i)
        for (std::size_t c = 0; c < Lchain.cols(); ++c)
            if (!std::isfinite(num_traits<T>::to_double(Lchain(i, c)))) Lchain(i, c) = zero;
    std::vector<int> Nchain(sn.nchains, 0);
    for (std::size_t c = 0; c < sn.nchains; ++c)
        Nchain[c] = static_cast<int>(std::llround(d.Nchain[c]));

    std::vector<std::size_t> infset;
    for (std::size_t i = 0; i < M; ++i)
        if (std::isinf(sn.stations[i].nservers)) infset.push_back(i);

    // G does not depend on how the delay stations are split: they aggregate by
    // the multinomial theorem, so it is taken with the infinite-server rows
    // summed into the think time.
    std::vector<bool> isinf(M, false);
    for (std::size_t k = 0; k < infset.size(); ++k) isinf[infset[k]] = true;
    std::size_t nq = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (!isinf[i]) ++nq;
    Matrix<T> Lq(nq, sn.nchains, zero);
    std::size_t a = 0;
    for (std::size_t i = 0; i < M; ++i) {
        if (isinf[i]) continue;
        for (std::size_t c = 0; c < sn.nchains; ++c) Lq(a, c) = Lchain(i, c);
        ++a;
    }
    Matrix<T> Z;
    if (!infset.empty()) {
        Z = Matrix<T>(1, sn.nchains, zero);
        for (std::size_t k = 0; k < infset.size(); ++k)
            for (std::size_t c = 0; c < sn.nchains; ++c) Z(0, c) += Lchain(infset[k], c);
    }
    const pfqn::NcResult<T> ca = pfqn::pfqn_ca(Lq, Nchain, Z);
    if (lG_out != nullptr) *lG_out = num_traits<T>::to_double(ca.lG);

    return pfqn::pfqn_jointmarg(nvec, Lchain, Nchain, infset, ca.G, engine,
                                static_cast<std::uint64_t>(opt.seed));
}

/**
 * Port of `@@SolverNC/getProbSysMarg.m`.
 *
 * Compare with `solver_nc_getprob_sys_aggr`, which fixes the PER-CLASS
 * population of every station and is a product form; each value returned here
 * is the sum of that one over every per-class table with these row sums.
 */
template <class T>
T solver_nc_getprob_sys_marg(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                             const std::vector<int>& nvec, const std::string& engine = "exact") {
    return solver_nc_jointmarg(sn, opt, nvec, engine, nullptr);
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_PROB_H
