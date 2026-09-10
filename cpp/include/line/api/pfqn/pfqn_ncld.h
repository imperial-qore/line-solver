/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_NCLD_H
#define LINE_API_PFQN_NCLD_H

/**
 * Normalizing constant of a LOAD-DEPENDENT closed network: the dispatcher.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ncld.m. The model reduction is the
 * same one pfqn_nc performs -- drop the empty classes, rescale per class, drop
 * the demand-free stations, peel off the classes confined to the delay -- with
 * one addition: the rate lattice mu follows the stations through the station
 * filter, so a dropped station takes its rates with it.
 *
 * DELAY FOLDING. pfqn_gld and pfqn_lldsingle take no think-time argument: an
 * infinite server is an ordinary row whose rate lattice is mu(i,k) = k, for
 * which the factorials cancel. When the model has a delay this port appends
 * one such row per think-time row, exactly as the reference does with
 * `Lz = [L;Z]; muz = [mu; repmat(1:size(mu,2),D,1)]`.
 *
 * DISPATCH. The exact ladder, which is the reference's `exact` branch and its
 * `default` branch whenever the Choudhury-Leung-Whitt gate declines:
 *
 *   R == 1                       -> pfqn_lldsingle
 *   M == 1 with a delay          -> pfqn_comomrm_ld
 *   M == 1 without a delay       -> pfqn_comomrm_ld with a zero think time
 *   otherwise                    -> pfqn_gld
 *
 * and, beside it, the ladder that answers with a LOGARITHM, mirroring pfqn_nc's
 * split between an exact `finish` and an estimator `finish_log`:
 *
 *   clw        -> pfqn_clw_lld    (generating-function inversion)
 *   is         -> pfqn_ld_is      (sample-an-ordering importance sampling)
 *   panald  -> pfqn_panaceald  (Mitra-McKenna asymptotic expansion)
 *   rd         -> pfqn_rd         (recursive decomposition)
 *   nrp / nrl  -> pfqn_nrp / pfqn_nrl (Norlund-Rice probit / logit)
 *   comomld    -> pfqn_comomrm_ld, or pfqn_rd where CoMoM-LD does not apply
 *   divdiff    -> pfqn_explicit_ld (divided difference over the LLD closed form)
 *
 * THE CLW GATE ON `default` IS A COST MODEL, NOT AN ACCURACY ONE. CLW and the
 * exact convolution compute the same constant, so the gate -- at most 5 classes,
 * at most 200 jobs, at most 2e7 predicted contour points -- only decides which
 * is cheaper. It is consulted in transcendental arithmetic only; in an exact
 * field `default` goes straight to the exact ladder, which changes the SPEED on
 * a subset of models and never the value.
 *
 * Arithmetic: the exact ladder and the whole reduction are EXACT-CAPABLE. The
 * log-domain ladder is refused by name outside transcendental arithmetic, on
 * pfqn_nc's grounds: a contour inversion, a Monte Carlo estimate and a
 * truncated asymptotic series have no meaning in the rational field, and
 * silently substituting the exact convolution would answer with an algorithm
 * the caller did not ask for.
 */

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_clw.h"
#include "line/api/pfqn/pfqn_comomrm_ld.h"
#include "line/api/pfqn/pfqn_explicit_ld.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_gldsingle.h"
#include "line/api/pfqn/pfqn_lldsingle.h"
#include "line/api/pfqn/pfqn_ld_is.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/api/pfqn/pfqn_nc.h"
#include "line/api/pfqn/pfqn_nre.h"
#include "line/api/pfqn/pfqn_nrl.h"
#include "line/api/pfqn/pfqn_panaceald.h"
#include "line/api/pfqn/pfqn_rd.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** The load-dependent methods this port dispatches. */
enum class NcldMethod {
    Default, Exact, Is, Clw, Panald, Rd, Nrp, Nrl, Nre, Comomld, Divdiff
};

inline const char* ncld_method_name(NcldMethod m) {
    switch (m) {
        case NcldMethod::Default: return "default";
        case NcldMethod::Exact: return "exact";
        case NcldMethod::Is: return "is";
        case NcldMethod::Clw: return "clw";
        case NcldMethod::Panald: return "panald";
        case NcldMethod::Rd: return "rd";
        case NcldMethod::Nrp: return "nrp";
        case NcldMethod::Nrl: return "nrl";
        case NcldMethod::Nre: return "nre";
        case NcldMethod::Comomld: return "comomld";
        case NcldMethod::Divdiff: return "divdiff";
    }
    return "default";
}

/** Map a method name to its enum; throws UnsupportedError on an unknown one. */
inline NcldMethod ncld_method_of(const std::string& s) {
    if (s == "default") return NcldMethod::Default;
    if (s == "exact") return NcldMethod::Exact;
    if (s == "is") return NcldMethod::Is;
    if (s == "clw") return NcldMethod::Clw;
    if (s == "pana" || s == "panald") return NcldMethod::Panald;
    if (s == "rd") return NcldMethod::Rd;
    if (s == "nrp") return NcldMethod::Nrp;
    if (s == "nrl") return NcldMethod::Nrl;
    if (s == "nre") return NcldMethod::Nre;
    if (s == "comomld") return NcldMethod::Comomld;
    if (s == "divdiff") return NcldMethod::Divdiff;
    throw UnsupportedError("pfqn_ncld: unrecognized method for solving load-dependent models: '" +
                           s + "'");
}

/** Refuse a load-dependent method in an arithmetic it has no meaning in. */
inline void pfqn_ncld_refuse(const std::string& method) {
    throw UnsupportedError(
        "pfqn_ncld: method '" + method +
        "' inverts a generating function, samples, or truncates an asymptotic series, and needs "
        "transcendental arithmetic. Use 'default', 'exact' or 'comomld' for an exact constant.");
}

template <class T>
struct NcldResult {
    T G;                 ///< normalizing constant
    double lG;           ///< its logarithm
    std::string method;  ///< the algorithm actually used
};

/**
 * @param L  (M x R) service demands
 * @param N  (R) populations, finite and nonnegative
 * @param Z  (K x R) think times
 * @param mu (M x >=Nt) load-dependent rate lattice
 * @param method requested algorithm
 * @param atol threshold below which a demand counts as zero; 0 for exact
 * @param nopt sample count, seed and tolerance the log-domain ladder reads
 */
template <class T>
NcldResult<T> pfqn_ncld(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                        const Matrix<T>& mu, NcldMethod method, const T& atol,
                        const NcOptions& nopt) {
    const std::size_t R0 = N.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    NcldResult<T> res;
    res.G = one;
    res.lG = 0.0;
    res.method = ncld_method_name(method);

    long Ntot = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_ncld: negative population");
        Ntot += v;
    }
    if (Ntot == 0) return res;
    const std::size_t M0 = L.empty() ? 0 : L.rows();
    if (M0 > 0 && L.cols() != R0)
        throw InputError("pfqn_ncld: L and N disagree on the class count");
    if (M0 > 0 && mu.rows() != M0)
        throw InputError("pfqn_ncld: mu and L disagree on the station count");
    if (M0 > 0 && static_cast<long>(mu.cols()) < Ntot)
        throw InputError("pfqn_ncld: mu has fewer rate columns than the total population");
    const std::size_t Nt = static_cast<std::size_t>(Ntot);

    // ---- drop the empty classes ----------------------------------------------
    std::vector<std::size_t> nnz;
    for (std::size_t r = 0; r < R0; ++r)
        if (N[r] > 0) nnz.push_back(r);
    const std::size_t R1 = nnz.size();

    // ---- rescale each class ---------------------------------------------------
    std::vector<T> scalevec(R1, one);
    Matrix<T> L1(M0, R1), Z1(Z.empty() ? 0 : Z.rows(), R1);
    for (std::size_t k = 0; k < R1; ++k) {
        const std::size_t r = nnz[k];
        T mx = zero;
        for (std::size_t i = 0; i < M0; ++i)
            if (L(i, r) > mx) mx = L(i, r);
        for (std::size_t i = 0; i < Z1.rows(); ++i)
            if (Z(i, r) > mx) mx = Z(i, r);
        if (mx > zero) scalevec[k] = mx;
        for (std::size_t i = 0; i < M0; ++i) L1(i, k) = L(i, r) / scalevec[k];
        for (std::size_t i = 0; i < Z1.rows(); ++i) Z1(i, k) = Z(i, r) / scalevec[k];
    }
    T Gscale = one;
    for (std::size_t k = 0; k < R1; ++k)
        Gscale *= num_pow_int(scalevec[k], static_cast<unsigned>(N[nnz[k]]));

    // ---- drop the stations with no demand, taking their rates with them -------
    std::vector<std::size_t> demSt;
    for (std::size_t i = 0; i < M0; ++i) {
        T rs = zero;
        for (std::size_t k = 0; k < R1; ++k) rs += L1(i, k);
        if (rs > atol) demSt.push_back(i);
    }
    const std::size_t M = demSt.size();
    Matrix<T> L2(M, R1), mu2(M, Nt);
    for (std::size_t a = 0; a < M; ++a) {
        for (std::size_t k = 0; k < R1; ++k) L2(a, k) = L1(demSt[a], k);
        for (std::size_t k = 0; k < Nt; ++k) mu2(a, k) = mu(demSt[a], k);
    }

    std::vector<int> N2(R1, 0);
    for (std::size_t k = 0; k < R1; ++k) N2[k] = N[nnz[k]];

    const auto delayG = [&](const std::vector<std::size_t>& cls) {
        T g = one;
        for (std::size_t k : cls) {
            T zs = zero;
            for (std::size_t i = 0; i < Z1.rows(); ++i) zs += Z1(i, k);
            g *= num_pow_int(zs, static_cast<unsigned>(N2[k])) /
                 num_factorial<T>(static_cast<unsigned>(N2[k]));
        }
        return g;
    };
    const auto finish = [&](const T& gcore) {
        res.G = Gscale * gcore;
        res.lG = num_traits<T>::log_as_double(res.G);
        return res;
    };
    // The log-domain finish, for pfqn_nc's reason: an estimator answers with lG,
    // and exponentiating it to take its logarithm again loses the answer once lG
    // passes ~709.
    const auto finish_log = [&](double lgcore) {
        res.lG = num_traits<T>::log_as_double(Gscale) + lgcore;
        res.G = num_traits<T>::from_double(std::exp(res.lG));
        return res;
    };

    T Ztot = zero;
    for (std::size_t i = 0; i < Z1.rows(); ++i)
        for (std::size_t k = 0; k < R1; ++k) Ztot += Z1(i, k);
    T Lsum = zero;
    for (std::size_t a = 0; a < M; ++a)
        for (std::size_t k = 0; k < R1; ++k) Lsum += L2(a, k);

    if (M == 0 || !(Lsum > atol)) {
        std::vector<std::size_t> all(R1);
        for (std::size_t k = 0; k < R1; ++k) all[k] = k;
        return finish(Ztot > atol ? delayG(all) : one);
    }
    if (M == 1 && !(Ztot > atol)) {
        // Single load-dependent station, no delay:
        // G = (sum N)! / prod N_r! * prod L_r^{N_r} / prod_{k<=Nt} mu(k).
        long tot = 0;
        for (int v : N2) tot += v;
        T g = num_factorial<T>(static_cast<unsigned>(tot));
        for (std::size_t k = 0; k < R1; ++k)
            g *= num_pow_int(L2(0, k), static_cast<unsigned>(N2[k])) /
                 num_factorial<T>(static_cast<unsigned>(N2[k]));
        for (std::size_t k = 0; k < Nt; ++k) {
            if (mu2(0, k) == zero) throw NumericError("pfqn_ncld: a load-dependent rate is zero");
            g /= mu2(0, k);
        }
        return finish(g);
    }

    // ---- classes confined to the delay ----------------------------------------
    std::vector<std::size_t> zdem, nzdem;
    for (std::size_t k = 0; k < R1; ++k) {
        T s = zero;
        for (std::size_t a = 0; a < M; ++a) s += L2(a, k);
        (s > atol ? nzdem : zdem).push_back(k);
    }
    const T Gzdem = zdem.empty() ? one : delayG(zdem);

    const std::size_t Rc = nzdem.size();
    Matrix<T> L3(M, Rc), Z3(Z1.rows(), Rc);
    std::vector<int> N3(Rc, 0);
    for (std::size_t a = 0; a < Rc; ++a) {
        for (std::size_t i = 0; i < M; ++i) L3(i, a) = L2(i, nzdem[a]);
        for (std::size_t i = 0; i < Z1.rows(); ++i) Z3(i, a) = Z1(i, nzdem[a]);
        N3[a] = N2[nzdem[a]];
    }
    T Z3tot = zero;
    for (std::size_t i = 0; i < Z3.rows(); ++i)
        for (std::size_t a = 0; a < Rc; ++a) Z3tot += Z3(i, a);

    // ---- the ladder that answers with a logarithm -----------------------------
    // The reference's cost model for 'default': CLW is preferred over the exact
    // convolution only where it is cheaper, so the three gates are class count,
    // total population and predicted contour points. `clw_defaults` supplies the
    // same inner lattice l_j the cost is predicted from, so the two cannot drift.
    long Nsum3 = 0;
    for (int v : N3) Nsum3 += v;
    bool default_clw = false;
    // NOT CONSULTED IN AN EXACT FIELD, and that is what keeps `default` working
    // there: the gate only decides which of two routes to the SAME constant is
    // cheaper, so an exact backend takes the exact ladder rather than being
    // refused for asking for a contour inversion it never named. Every caller of
    // the four-argument overload -- pfqn_ncldmx among them -- arrives on
    // `default`, so gating this at run time would refuse them under --arith
    // exact.
    if constexpr (num_traits<T>::has_transcendental) {
        if (method == NcldMethod::Default && M > 1 && Rc >= 2 && Rc <= 5 && Nsum3 <= 200) {
            std::vector<int> lat;
            std::vector<double> gam;
            // This cost estimate runs BEFORE any decomposition plan exists, so
            // it asks for the defaults of the identity ordering. clw_defaults
            // used to key the special cases on the chain POSITION and now keys
            // them on `depth`; depth[j] = j+1 is the same thing for an unsplit
            // ordering, so this reproduces the pre-decomposition defaults
            // (l = 1, 2, 2, 3, ...) exactly.
            std::vector<std::size_t> keep(Rc), depth(Rc);
            for (std::size_t a = 0; a < Rc; ++a) {
                keep[a] = a;
                depth[a] = a + 1;
            }
            detail::clw_defaults(Rc, keep, depth, ClwOptions(), lat, gam);
            double cost = 1.0;
            for (std::size_t a = 0; a < Rc; ++a)
                cost *= 2.0 * static_cast<double>(lat[a]) * static_cast<double>(N3[a]);
            default_clw = cost <= 2e7;  // ~2s at ~1e7 points/s
        }
    }
    // CoMoM-LD carries only the delay-plus-identical-stations shape; outside it
    // the reference warns and runs 'rd'. The substitution is visible in
    // res.method, which is what the reference's warning conveys.
    const bool comomld_falls_back =
        method == NcldMethod::Comomld && M > 1 && Z3tot > num_traits<T>::from_double(1e-14);
    const bool logdomain = default_clw || comomld_falls_back || method == NcldMethod::Is ||
                           method == NcldMethod::Clw || method == NcldMethod::Panald ||
                           method == NcldMethod::Rd || method == NcldMethod::Nrp ||
                           method == NcldMethod::Nrl || method == NcldMethod::Nre ||
                           method == NcldMethod::Divdiff;
    if (logdomain) {
        if constexpr (!num_traits<T>::has_transcendental) {
            // `default` never reaches here: the gate above is compiled out, so
            // every name arriving is one the caller asked for explicitly. A
            // `comomld` that fell back is named WITH its fallback, since the
            // caller's own word is exact-capable and the algorithm it resolved
            // to is not.
            pfqn_ncld_refuse(
                comomld_falls_back
                    ? std::string("comomld (which falls back to 'rd' on a multi-station model "
                                  "with a delay, CoMoM-LD carrying only the delay-plus-identical-"
                                  "stations shape)")
                    : std::string(ncld_method_name(method)));
        } else {
            const double lgz = num_traits<T>::log_as_double(Gzdem);
            // The aggregate think time per class, the reference's sum(Z,1).
            std::vector<T> Zv(Rc, zero);
            for (std::size_t i = 0; i < Z3.rows(); ++i)
                for (std::size_t a = 0; a < Rc; ++a) Zv[a] += Z3(i, a);

            if (default_clw || method == NcldMethod::Clw) {
                res.method = "clw";
                return finish_log(
                    lgz + num_traits<T>::to_double(pfqn_clw_lld(L3, N3, Zv, mu2, ClwOptions()).lG));
            }
            if (method == NcldMethod::Is) {
                McRng rng(static_cast<std::uint64_t>(nopt.seed));
                res.method = "is";
                return finish_log(lgz + pfqn_ld_is(L3, N3, Zv, mu2, nopt.samples, rng).lG);
            }
            if (method == NcldMethod::Panald) {
                const PanaceaLdResult<T> pa = pfqn_panaceald(L3, N3, Zv, mu2);
                if (!pa.normalUsage)
                    // The reference raises here rather than returning the NaN,
                    // because normal usage is the DOMAIN of the expansion and
                    // not a numerical failure it could retry out of.
                    throw UnsupportedError(
                        std::string("pfqn_ncld: the 'panald' asymptotic expansion does not "
                                    "apply to this model: ") +
                        (pa.reason ? pa.reason : "the expansion declined") +
                        ". Use 'exact', 'clw' or an approximate load-dependent method instead");
                res.method = "panald";
                return finish_log(lgz + num_traits<T>::to_double(pa.lG));
            }
            if (method == NcldMethod::Nrp) {
                std::vector<T> Nv(Rc, zero);
                for (std::size_t a = 0; a < Rc; ++a) Nv[a] = num_traits<T>::from_int(N3[a]);
                res.method = "nrp";
                return finish_log(lgz + num_traits<T>::to_double(pfqn_nrp(L3, Nv, Zv, mu2)));
            }
            if (method == NcldMethod::Nrl) {
                std::vector<T> Nv(Rc, zero);
                for (std::size_t a = 0; a < Rc; ++a) Nv[a] = num_traits<T>::from_int(N3[a]);
                res.method = "nrl";
                return finish_log(lgz + num_traits<T>::to_double(pfqn_nrl(L3, Nv, Zv, mu2)));
            }
            if (method == NcldMethod::Nre) {
                std::vector<T> Nv(Rc, zero);
                for (std::size_t a = 0; a < Rc; ++a) Nv[a] = num_traits<T>::from_int(N3[a]);
                res.method = "nre";
                return finish_log(lgz + num_traits<T>::to_double(pfqn_nre(L3, Nv, Zv, mu2)));
            }
            if (method == NcldMethod::Divdiff) {
                // Divided-difference closed form with the limited load-dependent
                // kernel of Casale-Harrison-Ong (Perform. Eval. 2021), Theorem 1.
                // A think time would have to enter g_sigma, whose closed form
                // covers queues only, so it is refused here as pfqn_nc refuses it
                // in the fixed-rate case. The delay-CONFINED classes already left
                // in Gzdem, so lgz still applies.
                if (Z3tot > zero)
                    throw UnsupportedError(
                        "pfqn_ncld: the 'divdiff' method requires a model without think time, "
                        "which needs the integral form of Corollary 3.4. Use 'exact' or "
                        "'default'");
                const ExplicitResult<T> ex = pfqn_explicit_ld(L3, N3, mu2);
                // WARNINGS BECOME FLAGS on this side, so the loss is acted on here
                // rather than printed: the reference warns past 15 digits and hands
                // the number back, which is only safe because the user sees the
                // warning. With no such channel, returning a constant double
                // precision cannot carry would be a silent wrong answer.
                if (!ex.valid || ex.lossDigits > 15)
                    throw NumericError(
                        "pfqn_ncld: the 'divdiff' closed form was exhausted by cancellation on "
                        "this model (" + std::to_string(ex.lossDigits) +
                        " decimal digits lost). Use 'exact', 'comomld' or 'rd', or merge the "
                        "near-tied scaled demands with a looser tolerance");
                res.method = "divdiff.ld/" + ex.method;
                return finish_log(lgz + num_traits<T>::to_double(ex.lG));
            }
            // 'rd', and the CoMoM-LD fallback onto it.
            res.method = "rd";
            return finish_log(lgz + pfqn_rd(L3, N3, Z3, mu2).lGN);
        }
    }

    // An explicit 'comomld' that did NOT fall back goes straight to CoMoM-LD,
    // whatever the station count: the reference calls it unconditionally on this
    // side of the gate, where the exact ladder below would prefer pfqn_gld.
    if (method == NcldMethod::Comomld) {
        res.method = "comomld";
        return finish(Gzdem * pfqn_comomrm_ld(L3, N3, Z3, mu2).G);
    }

    // ---- fold the delay rows into the demand matrix ----------------------------
    const std::size_t D = Z3tot > atol ? Z3.rows() : 0;
    Matrix<T> Lz(M + D, Rc), muz(M + D, Nt);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t a = 0; a < Rc; ++a) Lz(i, a) = L3(i, a);
        for (std::size_t k = 0; k < Nt; ++k) muz(i, k) = mu2(i, k);
    }
    for (std::size_t d = 0; d < D; ++d) {
        for (std::size_t a = 0; a < Rc; ++a) Lz(M + d, a) = Z3(d, a);
        for (std::size_t k = 0; k < Nt; ++k)
            muz(M + d, k) = num_traits<T>::from_int(static_cast<long>(k) + 1);
    }

    T gcore = one;
    if (Rc == 1) {
        int n1 = N3[0];
        gcore = pfqn_lldsingle(Lz, n1, muz).G;
        res.method = "exact/gld";
    } else if (M == 1 && Z3tot > atol) {
        gcore = pfqn_comomrm_ld(L3, N3, Z3, mu2).G;
        res.method = "exact/comomld";
    } else if (M == 1) {
        gcore = pfqn_comomrm_ld(L3, N3, Matrix<T>(), mu2).G;
        res.method = "exact/comomld";
    } else {
        gcore = pfqn_gld(Lz, N3, muz).G;
        res.method = "exact/gld";
    }

    return finish(Gzdem * gcore);
}

/** Overload with the reference's default sample count, seed and tolerance. */
template <class T>
NcldResult<T> pfqn_ncld(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                        const Matrix<T>& mu, NcldMethod method, const T& atol) {
    return pfqn_ncld(L, N, Z, mu, method, atol, NcOptions());
}

/** Overload with the exact (zero-tolerance) filters. */
template <class T>
NcldResult<T> pfqn_ncld(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                        const Matrix<T>& mu) {
    return pfqn_ncld(L, N, Z, mu, NcldMethod::Default, num_traits<T>::from_int(0), NcOptions());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_NCLD_H
