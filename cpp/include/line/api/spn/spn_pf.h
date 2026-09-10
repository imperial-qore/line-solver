#ifndef LINE_API_SPN_SPN_PF_H
#define LINE_API_SPN_SPN_PF_H

/**
 * @file spn_pf.h
 * @brief Product form of a stochastic Petri net: decide whether one exists and
 *        derive the per-level factors g_l that `mdd_rec` and `spn_metrics` take
 *        as input.
 *
 * THIS IS THE PART THE MDD-REC PAPER DECLARES OUT OF SCOPE (FGCS Sec. 3.2).
 * Every other function in api/spn receives the g_l already formed; this one
 * derives them from the net, which is what lets a solver reach them.
 *
 * THE THEORY, IN ONE PARAGRAPH. Write I(t), O(t) for the input and output
 * vectors of mode t and lambda_t for its rate constant. Henderson-Taylor and
 * Coleman-Henderson-Taylor show that a net whose firing rate has the form
 *
 *   r_t(m) = lambda_t psi(m - I(t)) / psi(m),     m >= I(t)
 *
 * has invariant measure pi(m) = psi(m) prod_l y_l^{m_l} whenever the positive
 * vector y satisfies COMPLEX BALANCE: reading the distinct vectors appearing as
 * some I(t) or O(t) as the COMPLEXES of the net, the flow into every complex
 * equals the flow out of it,
 *
 *   sum_{t : O(t)=v} lambda_t y^{I(t)} = ( sum_{t : I(t)=v} lambda_t ) y^v.
 *
 * Two choices of psi are realisable in LINE's own rate law, and they are the two
 * tested for here:
 *
 *   psi = 1              r_t = lambda_t, the rate of a SINGLE-SERVER mode.
 *                        pi(m) = prod_l y_l^{m_l}, so g_l(k) = y_l^k.
 *   psi = prod_l 1/m_l!  r_t = lambda_t prod_l m_l!/(m_l-I_l)!, MASS ACTION,
 *                        reached through a marking-dependent firing rate or,
 *                        for a mode drawing one token from one place, by
 *                        infinite-server semantics.
 *                        pi(m) = prod_l y_l^{m_l}/m_l!, so g_l(k) = y_l^k/k!.
 *
 * Which one holds is not guessed from the model API: the effective rate LINE
 * would use, lambda_t min(enabling degree, servers) g(m), is EVALUATED at every
 * reachable marking and compared against both laws. A net that matches neither
 * under one common psi is refused by name, never approximated.
 *
 * SOLVING FOR y. Complex balance reads A_lambda Psi(y) = 0 with A_lambda the
 * Laplacian of the weighted digraph on complexes and Psi(y)_v = y^v. That
 * Laplacian is the TRANSPOSED GENERATOR of a Markov chain that hops from complex
 * to complex at the rate of the mode joining them, so its kernel on one linkage
 * class is that chain's stationary distribution and `ctmc_solve` returns it --
 * strictly positive exactly when the class is strongly connected, which is weak
 * reversibility. With that positive vector kappa in hand y follows from the
 * LINEAR system in x = log y,
 *
 *   (v - v0) x = log kappa_v - log kappa_v0,   v, v0 in the same linkage class.
 *
 * Feinberg's Deficiency Zero Theorem says this system is consistent for every
 * choice of rate constants when the net is weakly reversible and its deficiency
 * c - l - s is zero, which is why those two numbers are reported; but
 * consistency is CHECKED rather than assumed, so a net of positive deficiency
 * whose particular rates still admit a complex-balanced point is accepted on the
 * evidence.
 *
 * THE GAUGE, AND WHY THE MINIMUM-NORM SOLUTION IS THE CANONICAL ONE. Complex
 * balance fixes y only up to y -> y .* exp(u) for any u orthogonal to the
 * stoichiometric subspace S. Such a shift multiplies pi(m) by exp(u'm), which is
 * CONSTANT on one compatibility class, so every reported measure is invariant
 * under it -- but the normalising constant G itself is not, it scales by that
 * constant. A gauge must therefore be FIXED, or the four codebases would report
 * four different G on the same net. The one fixed here is x in the row space of
 * the constraint matrix, i.e. the minimum-norm solution, reached in a form that
 * is unique whichever least-squares primitive a codebase carries: solve
 * (rows rows^T) w = rhs and set x = rows^T w. Any two solutions w of that system
 * give the SAME rows^T w, so the answer does not depend on how the rank-deficient
 * solve breaks its tie.
 *
 * ARITHMETIC. The logarithm and the exponential are unavoidable here -- y is the
 * exponential of a least-squares solution -- so the whole derivation needs a
 * transcendental field and is refused under exact arithmetic by name. What
 * consumes the g_l afterwards, `mdd_rec` and `spn_metrics`, stays rational.
 *
 * References:
 *   J. L. Coleman, W. Henderson, P. G. Taylor, "Product form equilibrium
 *   distributions and a convolution algorithm for stochastic Petri nets",
 *   Performance Evaluation 26(3), 1996.
 *   M. Feinberg, "Complex balancing in general kinetic systems", Arch. Rational
 *   Mech. Anal. 49, 1972.
 *   D. F. Anderson, G. Craciun, T. G. Kurtz, "Product-form stationary
 *   distributions for deficiency zero chemical reaction networks", Bull. Math.
 *   Biol. 72, 2010.
 *
 * @see mdd_rec, spn_metrics, spn_mdd, spn_conv
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/spn/spn_mdd.h"
#include "line/util/error.h"
#include "line/util/lstsq.h"
#include "line/util/matrix.h"

namespace line {
namespace spn {

/** Options of the product-form derivation. */
struct SpnPfOptions {
    /** Per-place-level token bound, passed to spn_mdd; empty infers it. */
    std::vector<double> bound;
    /** Relative tolerance of the rate-law and complex-balance checks. */
    double tol = 1e-9;
    bool verbose = false;
};

/** The product form, and the certificate that it is one. */
template <class T>
struct SpnPfResult {
    /** g[l][k] = g_l(k), ready for mdd_rec. */
    std::vector<std::vector<T>> g;
    /** Positive vector solving complex balance. */
    std::vector<T> y;
    /** "geometric" or "massaction", the psi that was found. */
    std::string kind;
    /** The distinct complexes, one row each. */
    std::vector<std::vector<double>> complexes;
    int deficiency = 0;
    std::size_t linkage = 0;
    std::size_t srank = 0;
    bool weakly_reversible = false;
    /** Relative complex-balance residual at y. */
    double residual = 0.0;
    SpnResult<T> spn;
};

namespace detail {

/** min(enabling degree, servers): the sets of tokens firing at once. */
template <class T>
double spn_pf_servers(const std::vector<double>& m, const SpnMode<T>& mde, std::size_t L) {
    double deg = std::numeric_limits<double>::infinity();
    for (std::size_t l = 0; l < L; ++l)
        if (mde.enab[l] > 0) deg = std::min(deg, std::floor(m[l] / mde.enab[l]));
    if (!std::isfinite(deg)) deg = 1.0;   // consumes nothing: always one set
    return std::min(deg, mde.srv);
}

/** prod_l m_l!/(m_l - I_l)!, the ordered ways to pick the input tokens. */
inline double spn_pf_massaction(const std::vector<double>& m, const std::vector<double>& enab,
                                std::size_t L) {
    double r = 1.0;
    for (std::size_t l = 0; l < L; ++l)
        for (int j = 0; j < static_cast<int>(enab[l]); ++j) r *= (m[l] - j);
    return r;
}

inline bool spn_pf_close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

/** Index of one complex, appended in first-seen order so the indices agree with
 *  the MATLAB, Java and python twins. */
inline std::size_t spn_pf_complex(const std::vector<double>& v,
                                  std::vector<std::vector<double>>& clist,
                                  std::map<std::string, std::size_t>& seen) {
    std::string key;
    for (std::size_t l = 0; l < v.size(); ++l) key += std::to_string(v[l]) + ",";
    const std::map<std::string, std::size_t>::const_iterator it = seen.find(key);
    if (it != seen.end()) return it->second;
    clist.push_back(v);
    seen[key] = clist.size() - 1;
    return clist.size() - 1;
}

inline std::vector<bool> spn_pf_reach(std::size_t v0, const std::vector<std::size_t>& from,
                                      const std::vector<std::size_t>& to, std::size_t c) {
    std::vector<bool> seen(c, false);
    seen[v0] = true;
    std::vector<std::size_t> stack(1, v0);
    while (!stack.empty()) {
        const std::size_t u = stack.back();
        stack.pop_back();
        for (std::size_t e = 0; e < from.size(); ++e)
            if (from[e] == u && !seen[to[e]]) {
                seen[to[e]] = true;
                stack.push_back(to[e]);
            }
    }
    return seen;
}

inline std::string spn_pf_wrtext(bool wr) {
    return wr ? "weakly reversible" : "not weakly reversible";
}

}  // namespace detail

/**
 * Derive the product form of a stochastic Petri net.
 *
 * @param sn      the network structure of a net holding Places and Transitions
 * @param options tolerances and bounds
 * @return the per-level factors and the certificate
 */
template <class T>
SpnPfResult<T> spn_pf(const qn::NetworkStruct<T>& sn, const SpnPfOptions& options = SpnPfOptions()) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "spn_pf: y is the exponential of a least-squares solution of the log-complex-balance "
            "equations, so the derivation needs a transcendental field and is unavailable under "
            "exact arithmetic. Supply the g_l directly to mdd_rec / spn_metrics, which stay "
            "rational");
    } else {
        SpnOptions mddopt;
        mddopt.descriptor = false;
        mddopt.bound = options.bound;
        SpnResult<T> spn = spn_mdd<T>(sn, mddopt);
        const SpnInfo<T>& info = spn.info;

        const std::size_t L = info.nplacelevels;
        const std::vector<SpnMode<T>>& md = info.modes;
        const std::size_t E = md.size();
        if (E == 0) throw InputError("spn_pf: the net has no timed mode");

        // ---- a queueing place holds an embedded server, not a token container
        for (std::size_t pp = 0; pp < info.places.size(); ++pp) {
            const std::size_t nd = info.places[pp];
            const std::size_t ist = sn.nodes[nd - 1].station;
            if (ist >= 1 && ist <= sn.stations.size() &&
                sn.stations[ist - 1].sched != lang::SchedStrategy::INF)
                throw UnsupportedError(
                    "spn_pf: place " + info.placenames[pp] +
                    " is a QUEUEING place: its embedded service is state that the marking does "
                    "not carry, so the net is not the token-container Petri net this product form "
                    "is written for");
        }

        // ---- the rate constants and the structural vectors
        std::vector<double> lambda(E, 0.0);
        std::vector<std::vector<double>> Iv(E), Ov(E);
        for (std::size_t e = 0; e < E; ++e) {
            lambda[e] = num_traits<T>::to_double(md[e].D1(0, 0));
            Iv[e] = md[e].enab;
            Ov[e] = md[e].fire;
            bool consumes = false;
            for (std::size_t l = 0; l < L; ++l) {
                if (std::isfinite(md[e].inhib[l]))
                    throw UnsupportedError(
                        "spn_pf: mode " + std::to_string(md[e].mode + 1) + " of node " +
                        std::to_string(md[e].trans) +
                        " has an inhibitor arc. An inhibitor zeroes the firing rate on markings "
                        "that still satisfy m >= I(t), so the rate is not lambda*psi(m-I)/psi(m) "
                        "on any psi and the net has no product form of this kind");
                if (md[e].enab[l] > 0) consumes = true;
            }
            if (md[e].srv != 1 && !consumes)
                throw InputError("spn_pf: mode " + std::to_string(md[e].mode + 1) + " of node " +
                                 std::to_string(md[e].trans) +
                                 " has several servers but consumes from no place, so its "
                                 "enabling degree is unbounded and its firing rate undefined");
            if (!(lambda[e] > 0))
                throw InputError("spn_pf: mode " + std::to_string(md[e].mode + 1) + " of node " +
                                 std::to_string(md[e].trans) +
                                 " has a non-positive firing rate");
        }

        // ---- which psi does LINE's own rate law follow on this net?
        const std::vector<std::vector<int>> states = info.diagram.enumerate();
        bool okgeo = true, okma = true;
        for (std::size_t s = 0; s < states.size(); ++s) {
            std::vector<double> m(L, 0.0);
            for (std::size_t l = 0; l < L; ++l) m[l] = states[s][l];
            for (std::size_t e = 0; e < E; ++e) {
                bool enabled = true;
                for (std::size_t l = 0; l < L && enabled; ++l)
                    if (m[l] < md[e].enab[l]) enabled = false;
                if (!enabled) continue;
                double actual = lambda[e] * detail::spn_pf_servers<T>(m, md[e], L);
                if (md[e].dep) {
                    // The multiplier is written against the per-NODE, class-summed
                    // marking, which is what the CTMC hands it.
                    std::vector<T> mk(info.nnodes, num_traits<T>::from_int(0));
                    for (std::size_t pp = 0; pp < info.places.size(); ++pp)
                        mk[info.places[pp] - 1] = num_traits<T>::from_double(m[pp]);
                    actual *= num_traits<T>::to_double(md[e].dep(mk));
                }
                const double ma = lambda[e] * detail::spn_pf_massaction(m, md[e].enab, L);
                okgeo = okgeo && detail::spn_pf_close(actual, lambda[e], options.tol);
                okma = okma && detail::spn_pf_close(actual, ma, options.tol);
                if (!okgeo && !okma)
                    throw UnsupportedError(
                        "spn_pf: mode " + std::to_string(md[e].mode + 1) + " of node " +
                        std::to_string(md[e].trans) + " fires at rate " + std::to_string(actual) +
                        " in a reachable marking, which is neither its rate constant "
                        "(single-server, psi = 1) nor its mass-action rate " +
                        std::to_string(ma) +
                        " (psi = prod 1/m!). LINE's rate law on this mode is "
                        "lambda*min(enabling degree, servers)*g(m), and no psi puts that in the "
                        "form lambda*psi(m-I)/psi(m)");
            }
        }
        const std::string kind = okgeo ? "geometric" : "massaction";

        // ---- complexes and the weighted digraph on them
        std::vector<std::vector<double>> C;
        std::map<std::string, std::size_t> seen;
        std::vector<std::size_t> src(E), dst(E);
        for (std::size_t e = 0; e < E; ++e) src[e] = detail::spn_pf_complex(Iv[e], C, seen);
        for (std::size_t e = 0; e < E; ++e) dst[e] = detail::spn_pf_complex(Ov[e], C, seen);
        const std::size_t c = C.size();

        // ---- Laplacian: A(j,i) is the rate of the arc i -> j
        Matrix<double> A(c, c);
        for (std::size_t i = 0; i < c; ++i)
            for (std::size_t j = 0; j < c; ++j) A(i, j) = 0.0;
        for (std::size_t e = 0; e < E; ++e) {
            if (src[e] == dst[e]) continue;      // a mode that moves nothing
            A(dst[e], src[e]) += lambda[e];
            A(src[e], src[e]) -= lambda[e];
        }

        // ---- linkage classes, weak reversibility, deficiency
        std::vector<int> lclass(c, -1);
        std::size_t nlink = 0;
        for (std::size_t v = 0; v < c; ++v) {
            if (lclass[v] >= 0) continue;
            std::vector<std::size_t> stack(1, v);
            lclass[v] = static_cast<int>(nlink);
            while (!stack.empty()) {
                const std::size_t u = stack.back();
                stack.pop_back();
                for (std::size_t e = 0; e < E; ++e) {
                    std::size_t w = c;
                    if (src[e] == u) w = dst[e];
                    else if (dst[e] == u) w = src[e];
                    if (w < c && lclass[w] < 0) {
                        lclass[w] = static_cast<int>(nlink);
                        stack.push_back(w);
                    }
                }
            }
            ++nlink;
        }
        bool wr = true;
        for (std::size_t b = 0; b < nlink && wr; ++b) {
            std::size_t v0 = c;
            for (std::size_t v = 0; v < c && v0 == c; ++v)
                if (lclass[v] == static_cast<int>(b)) v0 = v;
            const std::vector<bool> fwd = detail::spn_pf_reach(v0, src, dst, c);
            const std::vector<bool> bwd = detail::spn_pf_reach(v0, dst, src, c);
            for (std::size_t v = 0; v < c; ++v)
                if (lclass[v] == static_cast<int>(b) && (!fwd[v] || !bwd[v])) wr = false;
        }
        Matrix<double> netm(E, L);
        for (std::size_t e = 0; e < E; ++e)
            for (std::size_t l = 0; l < L; ++l) netm(e, l) = Ov[e][l] - Iv[e][l];
        Matrix<double> netmc = netm;
        const std::size_t srank = rref(netmc, line::detail::lstsq_tolerance(netm)).size();
        const int deficiency = static_cast<int>(c) - static_cast<int>(nlink) -
                               static_cast<int>(srank);

        // ---- kappa: the positive balance flow on each linkage class
        std::vector<double> kappa(c, 0.0);
        for (std::size_t b = 0; b < nlink; ++b) {
            std::vector<std::size_t> idx;
            for (std::size_t v = 0; v < c; ++v)
                if (lclass[v] == static_cast<int>(b)) idx.push_back(v);
            if (idx.size() == 1) {
                kappa[idx[0]] = 1.0;
                continue;
            }
            // The block is the transposed generator of the complex-hopping chain,
            // so its kernel is that chain's stationary law.
            Matrix<double> Qb(idx.size(), idx.size());
            for (std::size_t i = 0; i < idx.size(); ++i)
                for (std::size_t j = 0; j < idx.size(); ++j) Qb(i, j) = A(idx[j], idx[i]);
            const std::vector<double> pb = mc::ctmc_solve<double>(Qb);
            double mx = 0.0;
            for (std::size_t i = 0; i < idx.size(); ++i) mx = std::max(mx, pb[i]);
            for (std::size_t i = 0; i < idx.size(); ++i) {
                if (!(pb[i] > 0))
                    throw UnsupportedError(
                        "spn_pf: linkage class " + std::to_string(b + 1) + " of the complex graph "
                        "carries no flow through complex " + std::to_string(idx[i] + 1) +
                        ", so the net admits no positive complex-balanced point. A weakly "
                        "reversible net has a strictly positive balance flow on every linkage "
                        "class; this one is " + detail::spn_pf_wrtext(wr));
                kappa[idx[i]] = pb[i] / mx;
            }
        }

        // ---- x = log y from the linear system on each linkage class
        std::vector<std::vector<double>> rows;
        std::vector<double> rhs;
        for (std::size_t b = 0; b < nlink; ++b) {
            std::size_t v0 = c;
            for (std::size_t v = 0; v < c; ++v) {
                if (lclass[v] != static_cast<int>(b)) continue;
                if (v0 == c) {
                    v0 = v;
                    continue;
                }
                std::vector<double> row(L, 0.0);
                for (std::size_t l = 0; l < L; ++l) row[l] = C[v][l] - C[v0][l];
                rows.push_back(row);
                rhs.push_back(std::log(kappa[v]) - std::log(kappa[v0]));
            }
        }
        std::vector<double> x(L, 0.0);
        if (!rows.empty()) {
            // minimum norm through the row space; see the gauge note in the header
            const std::size_t nr = rows.size();
            Matrix<double> RRt(nr, nr);
            for (std::size_t i = 0; i < nr; ++i)
                for (std::size_t j = 0; j < nr; ++j) {
                    double s = 0.0;
                    for (std::size_t l = 0; l < L; ++l) s += rows[i][l] * rows[j][l];
                    RRt(i, j) = s;
                }
            const LstsqResult<double> w = lstsq<double>(RRt, rhs, line::detail::lstsq_tolerance(RRt));
            for (std::size_t l = 0; l < L; ++l) {
                double s = 0.0;
                for (std::size_t i = 0; i < nr; ++i) s += rows[i][l] * w.x[i];
                x[l] = s;
            }
            double res = 0.0, scale = 1.0;
            for (std::size_t i = 0; i < nr; ++i) {
                double s = 0.0;
                for (std::size_t l = 0; l < L; ++l) s += rows[i][l] * x[l];
                res = std::max(res, std::fabs(s - rhs[i]));
                scale = std::max(scale, std::fabs(rhs[i]));
            }
            if (res > options.tol * scale)
                throw UnsupportedError(
                    "spn_pf: the complex-balance equations are inconsistent (residual " +
                    std::to_string(res) + "): this net has no product form of the tested kind at "
                    "these rates. Its deficiency is " + std::to_string(deficiency) + " and it is " +
                    detail::spn_pf_wrtext(wr) + "; the Deficiency Zero Theorem guarantees a "
                    "solution only at deficiency 0 with weak reversibility");
        }
        std::vector<double> yd(L, 1.0);
        for (std::size_t l = 0; l < L; ++l) yd[l] = std::exp(x[l]);

        // ---- verify complex balance itself, which is what makes pi stationary
        std::vector<double> psi(c, 1.0);
        for (std::size_t v = 0; v < c; ++v) {
            double p = 1.0;
            for (std::size_t l = 0; l < L; ++l) p *= std::pow(yd[l], C[v][l]);
            psi[v] = p;
        }
        double resb = 0.0, scaleb = std::numeric_limits<double>::min();
        for (std::size_t v = 0; v < c; ++v) {
            double s = 0.0, sa = 0.0;
            for (std::size_t u = 0; u < c; ++u) {
                s += A(v, u) * psi[u];
                sa += std::fabs(A(v, u)) * psi[u];
            }
            resb = std::max(resb, std::fabs(s));
            scaleb = std::max(scaleb, sa);
        }
        if (resb > options.tol * scaleb)
            throw UnsupportedError(
                "spn_pf: complex balance fails at the computed point (relative residual " +
                std::to_string(resb / scaleb) + "), so the product form would not be stationary");

        // ---- the per-level factors, tabulated over the reachable domain
        SpnPfResult<T> out;
        out.g.resize(L);
        out.y.resize(L);
        for (std::size_t l = 0; l < L; ++l) {
            out.y[l] = num_traits<T>::from_double(yd[l]);
            const int d = spn.mdds.domain[l];
            out.g[l].assign(d, num_traits<T>::from_int(1));
            T fact = num_traits<T>::from_int(1), pw = num_traits<T>::from_int(1);
            for (int k = 0; k < d; ++k) {
                if (k > 0) {
                    fact = T(fact * num_traits<T>::from_int(k));
                    pw = T(pw * out.y[l]);
                }
                out.g[l][k] = kind == "massaction" ? T(pw / fact) : pw;
            }
        }
        out.kind = kind;
        out.complexes = C;
        out.deficiency = deficiency;
        out.linkage = nlink;
        out.srank = srank;
        out.weakly_reversible = wr;
        out.residual = resb / scaleb;
        out.spn = spn;
        return out;
    }
}

}  // namespace spn
}  // namespace line

#endif  // LINE_API_SPN_SPN_PF_H
