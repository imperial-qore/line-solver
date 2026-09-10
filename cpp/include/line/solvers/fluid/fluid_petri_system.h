/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_PETRI_SYSTEM_H
#define LINE_SOLVERS_FLUID_PETRI_SYSTEM_H

/**
 * The closures, the closed enabling term, the rate vector and the drift
 * Jacobian of a stochastic Petri net's fluid limit.
 *
 * Port of `matlab/src/solvers/FLD/fluid_petri_theta.m`,
 * `fluid_petri_rates.m` and `fluid_petri_jacobian.m`, cross-checked against
 * `jar/src/main/java/jline/solvers/fluid/petri/PetriSystem.java` and
 * `PetriClosures.java`.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

#include "line/solvers/fluid/fluid_petri_terms.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {
namespace petri {

/** The standard normal CDF, without a statistics dependency. */
inline double norm_cdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

/** The standard normal density. */
inline double norm_pdf(double z) {
    return std::exp(-0.5 * z * z) / std::sqrt(2.0 * 3.14159265358979323846);
}

/**
 * E[min(X,Y)] and dE/dE[X] for jointly normal X, Y.
 *
 * @return the pair (expectation, derivative with respect to E[X])
 */
inline std::pair<double, double> min_closure(double n, double c, double s2, double vc,
                                             double cov) {
    if (std::isinf(c)) return std::make_pair(n, 1.0);
    const double th2 = s2 - 2.0 * cov + vc;
    if (th2 <= 0.0) {
        // A degenerate pair: the min is the smaller of the two exactly. The band
        // matches the other three codebases, so two of them stopping either side
        // of the kink read the same indicator.
        if (c - n > petri_fine_tol() * std::max(1.0, std::fabs(n)))
            return std::make_pair(n, 1.0);
        return std::make_pair(c, 0.0);
    }
    const double th = std::sqrt(th2);
    const double al = (n - c) / th;
    const double phi = norm_pdf(al);
    const double p = 1.0 - norm_cdf(al);
    const double h = n * p + c * (1.0 - p) - th * phi;
    return std::make_pair(h, p);
}

/** The result of the many-argument closure. */
struct MinMulti {
    double h = 0.0;           ///< E[min(X_1,...,X_A,c)]
    std::vector<double> g;    ///< dH/dMU(a): the probability that arc a binds
    double v = 0.0;           ///< Var[min] before the cap
};

/**
 * Min-normal closure of E[min(X_1,...,X_A,c)] by Clark's (1961) recursion.
 *
 * The recursion is ORDER DEPENDENT, as Clark's approximation always is: only the
 * first two moments of the running min are kept. The order is the caller's, i.e.
 * increasing state coordinate, which the layout fixes.
 *
 * @param mu means of the arguments, already scaled by the arc weights
 * @param S  covariance of the arguments, symmetric positive semi-definite
 * @param c  deterministic cap (the mode's server count); infinite for none
 */
inline MinMulti minmulti_closure(const std::vector<double>& mu, const Matrix<double>& S,
                                 double c) {
    const std::size_t A = mu.size();
    MinMulti out;
    // A mode with no input arc is enabled at degree one, the convention the
    // exact engines use (solver_ssa_nrm's spnEnDegree, after_global_event).
    if (A == 0) {
        out.h = std::min(1.0, c);
        out.v = 0.0;
        return out;
    }
    std::vector<std::vector<double>> s(A, std::vector<double>(A, 0.0));
    if (S.rows() >= A && S.cols() >= A)
        for (std::size_t i = 0; i < A; ++i)
            for (std::size_t j = 0; j < A; ++j) s[i][j] = S(i, j);

    double mz = mu[0];
    double vz = s[0][0];
    std::vector<double> covz(A, 0.0);
    for (std::size_t i = 0; i < A; ++i) covz[i] = s[0][i];
    std::vector<double> p(A, 1.0);

    for (std::size_t k = 1; k < A; ++k) {
        double th2 = vz + s[k][k] - 2.0 * covz[k];
        // A covariance beyond the Cauchy-Schwarz bound is not admissible.
        if (th2 < 0.0) th2 = 0.0;
        const double th = std::sqrt(th2);
        double pk, mw, vw;
        std::vector<double> cw(A, 0.0);
        if (th > 0.0) {
            const double al = (mz - mu[k]) / th;
            const double Phi = norm_cdf(al);
            const double phi = norm_pdf(al);
            pk = 1.0 - Phi;
            mw = mz * pk + mu[k] * Phi - th * phi;
            const double e2 = (mz * mz + vz) * pk + (mu[k] * mu[k] + s[k][k]) * Phi -
                              (mz + mu[k]) * th * phi;
            vw = e2 - mw * mw;
            // Clark's moment match can leave a negative variance where the two
            // arguments are nearly identical; the min of two equal normals has
            // the variance of either, which the clamp restores.
            if (vw < 0.0) vw = 0.0;
            for (std::size_t i = 0; i < A; ++i) cw[i] = covz[i] * pk + s[i][k] * Phi;
        } else {
            pk = (mu[k] - mz > petri_fine_tol() * std::max(1.0, std::fabs(mz))) ? 1.0 : 0.0;
            mw = std::min(mz, mu[k]);
            vw = pk * vz + (1.0 - pk) * s[k][k];
            for (std::size_t i = 0; i < A; ++i) cw[i] = covz[i] * pk + s[i][k] * (1.0 - pk);
        }
        p[k] = pk;
        mz = mw;
        vz = vw;
        covz = cw;
    }

    out.v = vz;
    // The cap, by the two-argument closure itself: c is deterministic, so its
    // variance and its covariance with the running min are both zero.
    const std::pair<double, double> capped = min_closure(mz, c, vz, 0.0, 0.0);
    out.h = capped.first;
    out.g.assign(A, 0.0);
    double tail = capped.second;
    for (std::size_t a = A - 1; a >= 1; --a) {
        out.g[a] = (1.0 - p[a]) * tail;
        tail *= p[a];
    }
    out.g[0] = tail;
    return out;
}

/** The closed enabling term of every mode and its derivative. */
struct PetriTheta {
    std::vector<double> theta;
    std::vector<std::vector<std::size_t>> dslot;
    std::vector<std::vector<double>> dval;
    std::vector<double> dep;
    std::vector<std::vector<std::size_t>> depslot;
    std::vector<std::vector<double>> depval;
};

namespace system_detail {

/** Sigma(a,b) as the closure reads it: zero where the drift never asks. */
inline double sig(const PetriTerms& t, const std::vector<double>& s2, std::size_t a,
                  std::size_t b) {
    if (a >= t.pair_index.size() || b >= t.pair_index.size()) return 0.0;
    const std::ptrdiff_t i = t.pair_index[a][b];
    return (i < 0 || static_cast<std::size_t>(i) >= s2.size()) ? 0.0
                                                               : s2[static_cast<std::size_t>(i)];
}

/**
 * The marking-dependent firing multiplier and its gradient.
 *
 * g is a user function of the (nnodes x nclasses) marking, so it is evaluated at
 * the MEAN marking -- a first-order closure of g, the same order at which
 * SolverCTMC evaluates it per state and the only one available without the
 * distribution of the marking. Its gradient has no analytic form, so it is taken
 * by central differences.
 */
inline void dep_of(const PetriTerms& t, const PetriMode& md, const std::vector<double>& x,
                   PetriTheta& th, std::size_t j) {
    Matrix<double> mm(t.I, t.K, 0.0);
    for (std::size_t s = 0; s < t.nm; ++s) mm(t.coord_node[s], t.coord_class[s]) = x[s];
    th.dep[j] = md.dep(mm);
    for (std::size_t s = 0; s < t.nm; ++s) {
        const double h = std::max(1e-6 * std::fabs(x[s]), 1e-6);
        const std::size_t i = t.coord_node[s], k = t.coord_class[s];
        Matrix<double> mp = mm, mn = mm;
        mp(i, k) = mm(i, k) + h;
        mn(i, k) = std::max(0.0, mm(i, k) - h);
        const double hh = mp(i, k) - mn(i, k);
        if (hh <= 0) continue;
        const double d = (md.dep(mp) - md.dep(mn)) / hh;
        if (d != 0) {
            th.depslot[j].push_back(s);
            th.depval[j].push_back(d);
        }
    }
}

}  // namespace system_detail

/**
 * The closed enabling term of every mode, and its derivative.
 *
 *   theta_j = ( prod_b Phi((thr_b - m_b)/sd_b) ) * E[ min_a(m_a/w_a), c_j ]
 *
 * Both factors collapse to their first-order form at zero variance -- Phi becomes
 * the hard indicator and the min closure becomes min() -- so the mean-field limit
 * is one code path, not two.
 *
 * THE INHIBITOR GATE IS WHY A PETRI NET NEEDS A SMOOTHED CLOSURE AT ALL, quite
 * apart from accuracy: the indicator is a step, and a Newton solver has no
 * derivative to descend on a step.
 *
 * THE VARIANCES ARE UNKNOWNS, NOT FUNCTIONS OF X: s2 is pinned by its own
 * consistency row in the DAE, so the derivative is with respect to the MEANS
 * only.
 */
inline PetriTheta petri_theta(const PetriTerms& t, const std::vector<double>& x,
                              const std::vector<double>& s2_in) {
    const std::size_t nmod = t.modes.size();
    PetriTheta th;
    th.theta.assign(nmod, 0.0);
    th.dslot.assign(nmod, std::vector<std::size_t>());
    th.dval.assign(nmod, std::vector<double>());
    th.dep.assign(nmod, 1.0);
    th.depslot.assign(nmod, std::vector<std::size_t>());
    th.depval.assign(nmod, std::vector<double>());
    const std::vector<double> zeros(std::max<std::size_t>(t.npair, 1), 0.0);
    const std::vector<double>& s2 = s2_in.empty() ? zeros : s2_in;

    for (std::size_t j = 0; j < nmod; ++j) {
        const PetriMode& md = t.modes[j];
        const std::size_t A = md.arc_slot.size();
        std::vector<double> mu(A, 0.0);
        for (std::size_t a = 0; a < A; ++a) mu[a] = x[md.arc_slot[a]] / md.arc_w[a];
        Matrix<double> Sarg;
        if (md.closable && A > 0) {
            Sarg = Matrix<double>(A, A, 0.0);
            for (std::size_t a = 0; a < A; ++a)
                for (std::size_t b = 0; b < A; ++b)
                    Sarg(a, b) = system_detail::sig(t, s2, md.arc_slot[a], md.arc_slot[b]) /
                                 (md.arc_w[a] * md.arc_w[b]);
        }
        const MinMulti mm = minmulti_closure(mu, Sarg, md.c);

        const std::size_t nb = md.inh_slot.size();
        std::vector<double> gate(nb, 0.0), dgate(nb, 0.0);
        double ginh = 1.0;
        for (std::size_t b = 0; b < nb; ++b) {
            const double mb = x[md.inh_slot[b]];
            const double thr = md.inh_thr[b];
            const double vb = system_detail::sig(t, s2, md.inh_slot[b], md.inh_slot[b]);
            if (vb > 0) {
                const double sd = std::sqrt(vb);
                const double zb = (thr - mb) / sd;
                gate[b] = norm_cdf(zb);
                dgate[b] = -norm_pdf(zb) / sd;
            } else {
                gate[b] = (mb < thr - petri_fine_tol() * std::max(1.0, thr)) ? 1.0 : 0.0;
                dgate[b] = 0.0;
            }
            ginh *= gate[b];
        }

        // An arc and an inhibitor arc may share a coordinate, so accumulate.
        std::map<std::size_t, double> acc;
        for (std::size_t a = 0; a < A; ++a)
            acc[md.arc_slot[a]] += ginh * mm.g[a] / md.arc_w[a];
        for (std::size_t b = 0; b < nb; ++b) {
            double others;
            if (gate[b] != 0) {
                others = ginh / gate[b];
            } else {
                others = 1.0;
                for (std::size_t q = 0; q < nb; ++q)
                    if (q != b) others *= gate[q];
            }
            acc[md.inh_slot[b]] += mm.h * others * dgate[b];
        }
        for (std::map<std::size_t, double>::const_iterator it = acc.begin(); it != acc.end();
             ++it) {
            th.dslot[j].push_back(it->first);
            th.dval[j].push_back(it->second);
        }
        th.theta[j] = ginh * mm.h;

        if (md.dep) system_detail::dep_of(t, md, x, th, j);
    }
    return th;
}

/**
 * The rate of every event column.
 *
 *   kind 1  firing of mode j    single phase: rateBase*theta*dep
 *                               multi  phase: rateBase*y(j,h)
 *   kind 2  internal phase change             rateBase*y(j,h)
 *   kind 3  exogenous arrival                 a constant
 *   kind 4  firing of an IMMEDIATE mode       phi_j, an algebraic unknown
 *   kind 5  the server latch                  mu_j, a free-sign unknown
 */
inline std::vector<double> petri_rates(const PetriTerms& t, const std::vector<double>& x,
                                       const std::vector<double>& phi,
                                       const std::vector<double>& mu, const PetriTheta& th) {
    const std::size_t nmod = t.modes.size();
    std::vector<std::ptrdiff_t> imm_pos(nmod, -1), latch_pos(nmod, -1);
    for (std::size_t q = 0; q < t.imm_idx.size(); ++q)
        imm_pos[t.imm_idx[q]] = static_cast<std::ptrdiff_t>(q);
    for (std::size_t q = 0; q < t.latch_mode.size(); ++q)
        latch_pos[t.latch_mode[q]] = static_cast<std::ptrdiff_t>(q);

    std::vector<double> r(t.nev, 0.0);
    for (std::size_t e = 0; e < t.nev; ++e) {
        const int k = t.ev_kind[e];
        if (k == 3) {
            r[e] = t.rate_base[e];
        } else if (k == 4) {
            const std::ptrdiff_t at = imm_pos[static_cast<std::size_t>(t.ev_mode[e])];
            r[e] = (phi.empty() || at < 0) ? 0.0 : phi[static_cast<std::size_t>(at)];
        } else if (k == 5) {
            const std::ptrdiff_t at = latch_pos[static_cast<std::size_t>(t.ev_mode[e])];
            r[e] = (mu.empty() || at < 0) ? 0.0 : mu[static_cast<std::size_t>(at)];
        } else {
            const std::size_t j = static_cast<std::size_t>(t.ev_mode[e]);
            const PetriMode& md = t.modes[j];
            if (md.nph == 1)
                r[e] = t.rate_base[e] * th.theta[j] * th.dep[j];
            else
                r[e] = t.rate_base[e] * x[md.zblk[static_cast<std::size_t>(t.ev_phase[e])]];
        }
    }
    return r;
}

/**
 * Drift Jacobian A = D * dR/dX.
 *
 * This is what the Lyapunov equation of the linear noise approximation is
 * written about, so it has to be the derivative of the SAME rate vector
 * `petri_rates` returns: a covariance solved about an inconsistent Jacobian is
 * not the covariance of anything. THE VARIANCES ARE HELD.
 */
inline Matrix<double> petri_jacobian(const PetriTerms& t, const PetriTheta& th) {
    Matrix<double> Jr(std::max<std::size_t>(t.nev, 1), t.nstate, 0.0);
    for (std::size_t e = 0; e < t.nev; ++e) {
        const int k = t.ev_kind[e];
        if (k == 3 || k == 4 || k == 5) continue;
        const std::size_t j = static_cast<std::size_t>(t.ev_mode[e]);
        const PetriMode& md = t.modes[j];
        const double base = t.rate_base[e];
        if (md.nph == 1) {
            for (std::size_t q = 0; q < th.dslot[j].size(); ++q)
                Jr(e, th.dslot[j][q]) += base * th.dep[j] * th.dval[j][q];
            for (std::size_t q = 0; q < th.depslot[j].size(); ++q)
                Jr(e, th.depslot[j][q]) += base * th.theta[j] * th.depval[j][q];
        } else {
            const std::size_t zc = md.zblk[static_cast<std::size_t>(t.ev_phase[e])];
            Jr(e, zc) += base;
        }
    }
    return matmul(t.D, Jr);
}

}  // namespace petri
}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_PETRI_SYSTEM_H
