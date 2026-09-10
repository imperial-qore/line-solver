/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_CLOSURES_H
#define LINE_SOLVERS_FLUID_FLUID_CLOSURES_H

/**
 * The moment closures the fluid drift is built from: `fluid_min_closure.m`,
 * `fluid_capacity_closure.m`, `fluid_lld_scaling.m`, `fluid_share_closure.m` and
 * `fluid_gps_share.m`.
 *
 * WHY THESE ARE ONE HEADER AND NOT PART OF THE DRIFT. Each answers the same
 * question about a different non-linear term of the rate function: what is
 * E[f(X)] when X is not known but its first two moments are. The first-order
 * fluid closure answers f(E[X]) for all of them, which is why they collapse to
 * nothing when the variance is zero and why every routine below returns the
 * first-order value on that path rather than special-casing it at the caller.
 * The drift needs the VALUE, the covariance equation needs the DERIVATIVE, and
 * they must be the derivative OF that value or the Lyapunov solve is linearizing
 * a different drift than the one integrated; each routine returns both together
 * for exactly that reason.
 *
 * THREE DISTINCT NON-LINEARITIES, and they are not interchangeable:
 *   min(n, c)      the server capacity, closed by a normal marginal
 *                  (`fluid_min_closure`, and `fluid_capacity_closure` once a
 *                  load-dependent alpha(n) multiplies it)
 *   w_j x_j / sum  the capacity share of PS and DPS, a RATIO, closed by the
 *                  delta method (`fluid_share_closure`)
 *   1{x_k >= 1}    the backlog indicator of GPS, closed by ENUMERATING the 2^K
 *                  backlog patterns (`fluid_gps_share`)
 * The third is not a correction but the whole mechanism: with continuous mass
 * every class is always backlogged, so a first-order closure prices GPS at its
 * heavy-traffic constant regardless of load.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** A closure's value and its first two derivatives with respect to the first mean. */
struct ClosureValue {
    double h = 0.0;
    double dh = 0.0;
    double d2h = 0.0;
};

/** The standard normal cdf, without a statistics library. */
inline double closure_normcdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

/** The standard normal pdf. */
inline double closure_normpdf(double z) {
    return std::exp(-0.5 * z * z) / std::sqrt(2.0 * 3.14159265358979323846);
}

/**
 * Port of `fluid_min_closure.m`: E[min(X,Y)] for jointly normal X, Y, and its
 * derivative with respect to E[X].
 *
 * The published closure is Guenther, Stefanek and Bradley (EPEW/UKPEW 2012,
 * LNCS 7587:32-47, eq. 4) in its general two-population form, of which SolverFLD
 * uses only the specialisation Y = c with c the deterministic server count. The
 * general arguments are kept so this IS the published closure rather than one
 * instance of it.
 *
 * With theta = 0 the expression collapses to min(n,c) and the derivative to
 * 1{n < c}, so the first-order closure shares this code path exactly. The
 * derivative AT the kink is 0, the right derivative of min(), which is the
 * convention the strict inequality of the first-order branch already implies.
 */
inline ClosureValue fluid_min_closure(double n, double c, double s2, double vc = 0.0,
                                      double cov_nc = 0.0) {
    ClosureValue r;
    double th2 = s2 - 2.0 * cov_nc + vc;
    if (th2 < 0.0) th2 = 0.0;  // a covariance beyond Cauchy-Schwarz is not admissible
    if (!(th2 > 0.0) || std::isinf(c)) {
        // THE INDICATOR CARRIES A BAND, and it is a cross-codebase requirement rather
        // than a modelling choice. A saturated fluid fixed point sits exactly AT n = c,
        // and each engine's ODE stops on its own residual: MATLAB lands at 1.0004 and
        // this port at 1 - 1.8e-13 on the same model, so a strict `n < c` reads
        // saturated in one and unsaturated in the other. That flips a whole Jacobian row
        // between zero and unit, and with it the hyperbolicity verdict that decides
        // whether SolverFluid answers with 'minnormal' or falls back to the first-order
        // method. A population within FineTol of the server count IS at the kink.
        r.h = std::min(n, c);
        r.dh = ((c - n) > lang::GlobalConstants::FineTol * std::max(1.0, n)) ? 1.0 : 0.0;
        return r;
    }
    const double th = std::sqrt(th2);
    const double d = (n - c) / th;
    const double Phid = closure_normcdf(d);
    const double phid = closure_normpdf(d);
    r.h = n * (1.0 - Phid) + c * Phid - th * phid;
    r.dh = 1.0 - Phid;
    // min() is piecewise linear, so its second derivative is carried entirely by the
    // atom at X = Y; smoothing over the normal marginal turns that atom into the
    // density. It stays zero on the degenerate branch above, which smooths nothing.
    r.d2h = -phid / th;
    // A POPULATION IS NONNEGATIVE AND THE NORMAL MARGINAL IS NOT. For X >= 0
    // pathwise min(X,c) >= 0, and min() being concave Jensen puts
    // E[min(X,c)] <= min(E[X],c), so the value belongs to [0, min(n,c)]. The normal
    // marginal has no such support, and the mass it places below zero drags the
    // expectation out of that range once the mean falls to about one standard
    // deviation: n = 0, c = 1, th = 0.664 returns -0.019, which is a station that
    // CREATES work. Project onto the admissible range and take the derivatives of
    // the bound that binds, so the Jacobian still matches h.
    const double hi = std::min(n, c);
    const bool at_lo = r.h < 0.0;
    const bool at_hi = r.h > hi;
    r.h = std::min(std::max(r.h, 0.0), hi);
    if (at_lo) {
        r.dh = 0.0;
    } else if (at_hi) {
        r.dh = (n < c) ? 1.0 : 0.0;
    }
    if (at_lo || at_hi) r.d2h = 0.0;
    return r;
}

/**
 * Port of `fluid_lld_scaling.m`: the limited load-dependent multiplier alpha at a
 * CONTINUOUS population, by linear interpolation of the integer table.
 *
 * `sn.lldscaling(i,:)` is tabulated at populations 1..lldlimit and the discrete
 * solvers read it at `min(n, lldlimit)`. The fluid state is continuous, so the
 * table is interpolated between consecutive entries and clamped to the first
 * entry below n = 1 and to the last above the table end, which is the clamping
 * the CTMC already applies.
 */
inline ClosureValue fluid_lld_scaling(const std::vector<double>& lldrow, double n) {
    ClosureValue r;
    if (lldrow.empty()) {
        r.h = 1.0;
        r.dh = 0.0;
        return r;
    }
    const std::size_t L = lldrow.size();
    if (n <= 1.0) {
        r.h = lldrow[0];
        return r;
    }
    if (n >= static_cast<double>(L)) {
        r.h = lldrow[L - 1];
        return r;
    }
    const std::size_t k = static_cast<std::size_t>(std::floor(n));  // 1 <= k <= L-1
    const double frac = n - static_cast<double>(k);
    r.h = lldrow[k - 1] * (1.0 - frac) + lldrow[k] * frac;
    r.dh = lldrow[k] - lldrow[k - 1];
    return r;
}

namespace detail {

/** psi(u) = A + B*u + C*u^2 on [p,q], from base(u)*alpha(u). */
struct PsiSegment {
    double A = 0.0, B = 0.0, C = 0.0;
};

/**
 * `local_segment` of `fluid_capacity_closure.m`.
 *
 * alpha is constant on [0,1] and on [L,inf) and linear on each unit interval, so
 * psi = base*alpha is piecewise quadratic with breakpoints at the integers and at
 * the saturation point; the breakpoint list guarantees the segment lies entirely
 * on one side of c, which is what lets base be either u or c and not both.
 */
inline PsiSegment psi_segment(double p, double q, double c, const std::vector<double>& lldrow,
                              bool is_inf) {
    const std::size_t L = lldrow.size();
    double a0 = 0.0, a1 = 0.0;
    if (p >= static_cast<double>(L)) {
        a0 = lldrow[L - 1];
    } else if (p < 1.0) {
        a0 = lldrow[0];
    } else {
        const std::size_t k = static_cast<std::size_t>(std::floor(p));
        a1 = lldrow[k] - lldrow[k - 1];
        a0 = lldrow[k - 1] - a1 * static_cast<double>(k);
    }
    PsiSegment s;
    if (is_inf || !std::isfinite(c) || q <= c) {
        s.A = 0.0;
        s.B = a0;
        s.C = a1;  // base = u
    } else {
        s.A = c * a0;
        s.B = c * a1;
        s.C = 0.0;  // base = c
    }
    return s;
}

/** Truncated moments E[X^j 1{p < X < q}] for X ~ Normal(n, s^2), j = 0,1,2. */
struct TruncMoments {
    double M0 = 0.0, M1 = 0.0, M2 = 0.0;
};

inline TruncMoments trunc_moments(double p, double q, double n, double s) {
    const double zp = (p - n) / s;
    const double Pp = closure_normcdf(zp);
    const double pp = closure_normpdf(zp);
    const bool qinf = std::isinf(q);
    const double Pq = qinf ? 1.0 : closure_normcdf((q - n) / s);
    const double pq = qinf ? 0.0 : closure_normpdf((q - n) / s);
    TruncMoments m;
    m.M0 = Pq - Pp;
    m.M1 = n * m.M0 + s * (pp - pq);
    m.M2 = qinf ? (n * n + s * s) * m.M0 + s * ((p + n) * pp)
                : (n * n + s * s) * m.M0 + s * ((p + n) * pp - (q + n) * pq);
    return m;
}

/** psi and its derivative at a continuous population, extended by zero below 0. */
inline ClosureValue psi_point(double u, double c, const std::vector<double>& lldrow, bool is_inf) {
    const ClosureValue a = fluid_lld_scaling(lldrow, u);
    const double base = is_inf ? u : std::min(u, c);
    const double dbase = is_inf ? 1.0 : ((u < c) ? 1.0 : 0.0);
    ClosureValue r;
    r.h = base * a.h;
    r.dh = dbase * a.h + base * a.dh;
    // base and alpha are both piecewise linear, so psi is piecewise quadratic and
    // psi'' = 2*base'*alpha' away from the breakpoints. The atoms AT the
    // breakpoints are not representable without a marginal to smooth them, and the
    // first-order closure this branch serves does not smooth them either.
    r.d2h = 2.0 * dbase * a.dh;
    if (u <= 0.0) {
        r.h = 0.0;
        r.dh = 0.0;
        r.d2h = 0.0;
    }
    return r;
}

}  // namespace detail

/**
 * Port of `fluid_capacity_closure.m`: E[psi(X)] and its derivative, where
 * psi(n) = min(n,c)*alpha(n) at a queueing station and n*alpha(n) at an infinite
 * server.
 *
 * THE LOAD-DEPENDENT INTEGRATION IS EXACT AND NOT A QUADRATURE, which matters
 * beyond accuracy: Gauss-Hermite on a fixed node set does not smooth the kinks of
 * a piecewise-linear integrand, it RELOCATES them, and the resulting E[psi] is
 * itself piecewise linear in n. Its second derivative is then zero almost
 * everywhere and `fluid_refine_meanfield` silently returns a null correction.
 * The segment-wise closed form below is what avoids that.
 *
 * The zero floor on the no-load-dependence branch is the reference's minimal
 * repair for the normal marginal putting mass below zero, where min(X,c) = X < 0:
 * once n is small against the standard deviation E[min(X,c)] itself goes
 * negative, which is a negative service rate and mass destroyed by the
 * integrator's non-negativity clamp.
 */
inline ClosureValue fluid_capacity_closure(double n, double c, double s2,
                                           const std::vector<double>& lldrow, bool is_inf) {
    if (lldrow.empty()) {
        ClosureValue r;
        if (is_inf) {
            r.h = n;
            r.dh = 1.0;
            return r;
        }
        r = fluid_min_closure(n, c, s2);
        if (r.h < 0.0) {
            r.h = 0.0;
            r.dh = 0.0;
            r.d2h = 0.0;
        }
        return r;
    }

    const std::size_t L = lldrow.size();
    if (!(s2 > 0.0)) return detail::psi_point(n, c, lldrow, is_inf);

    const double s = std::sqrt(s2);
    // The breakpoints of psi: the lattice of the table, the origin, and the
    // saturation point c when it falls beyond the table.
    std::vector<double> bps;
    for (std::size_t k = 0; k <= L; ++k) bps.push_back(static_cast<double>(k));
    if (!is_inf && std::isfinite(c) && c > static_cast<double>(L)) bps.push_back(c);
    std::sort(bps.begin(), bps.end());
    bps.erase(std::unique(bps.begin(), bps.end()), bps.end());

    ClosureValue r;
    // psi is extended by zero below the first breakpoint, so psi' jumps there too
    // and that atom belongs in psi'' exactly like the interior ones
    double Bprev = 0.0, Cprev = 0.0;
    for (std::size_t k = 0; k < bps.size(); ++k) {
        const double p = bps[k];
        const double q = (k + 1 < bps.size()) ? bps[k + 1] : std::numeric_limits<double>::infinity();
        const detail::PsiSegment g = detail::psi_segment(p, q, c, lldrow, is_inf);
        const detail::TruncMoments mo = detail::trunc_moments(p, q, n, s);
        r.h += g.A * mo.M0 + g.B * mo.M1 + g.C * mo.M2;
        r.dh += g.B * mo.M0 + 2.0 * g.C * mo.M1;
        r.d2h += 2.0 * g.C * mo.M0;
        // psi' jumps across this breakpoint, so psi'' carries an atom there; the
        // segment sum above sees only the quadratic part and would miss it
        const double jump = (g.B + 2.0 * g.C * p) - (Bprev + 2.0 * Cprev * p);
        r.d2h += jump * closure_normpdf((p - n) / s) / s;
        Bprev = g.B;
        Cprev = g.C;
    }
    return r;
}

/**
 * Port of `local_project_rate` in `ode_rates_closing_factors.m`: project a jointly
 * closed per-coordinate service share onto the set it has to live in, namely
 * `r >= 0`, `r <= xb` where that bound applies, and `sum(r) = tot`.
 *
 * THE JOINT CLOSURE IS AN EXPANSION AND CAN LEAVE THAT SET. `r = s*psi +
 * psi'*Cov(S,N)` adds a term that sums to ZERO over the coordinates, so it moves
 * mass between them and its entries can push one past either bound; the
 * first-order share `x_j/n_i*psi` cannot, being x_j scaled by `psi/n_i <= 1`.
 * Either breach ends the same way, because the integrator holds every coordinate
 * non-negative: `r_j > x_j` drains coordinate j faster than it holds, the state
 * goes negative and the clamp INJECTS mass.
 *
 * THE UPPER BOUND HOLDS ONLY WITHOUT LOAD DEPENDENCE, which is what CAPPED
 * selects. r is an expected NUMBER in service so `r_j <= x_j`, but
 * `psi(n) = min(n,c)*alpha(n)` folds the load-dependent scaling into the same
 * variable, and with alpha > 1 the first-order share itself exceeds x_j.
 *
 * Clip, then move the residual onto the coordinates that still have slack in
 * proportion to it, so `sum(r) = tot` survives and the station still clears what
 * its capacity closure says it clears. It is a NO-OP whenever the expansion stayed
 * inside the set, which is why models already inside it are bit-identical.
 */
inline void fluid_project_rate(std::vector<double>& r, const std::vector<double>& xb, bool capped,
                               double tot) {
    const double zt = lang::GlobalConstants::Zero;
    bool inside = true;
    for (std::size_t j = 0; j < r.size() && inside; ++j) {
        if (!(r[j] >= -zt)) inside = false;
        if (capped && !(r[j] <= xb[j] + zt)) inside = false;
    }
    if (inside) return;
    for (std::size_t j = 0; j < r.size(); ++j) {
        r[j] = std::max(r[j], 0.0);
        if (capped) r[j] = std::min(r[j], xb[j]);
    }
    for (std::size_t it = 0; it <= r.size(); ++it) {
        double sum = 0.0;
        for (std::size_t j = 0; j < r.size(); ++j) sum += r[j];
        const double d = tot - sum;
        if (std::fabs(d) <= zt) break;
        std::vector<double> slack(r.size(), 0.0);
        for (std::size_t j = 0; j < r.size(); ++j)
            slack[j] = (d > 0.0) ? (capped ? (xb[j] - r[j]) : 1.0) : r[j];
        double tsl = 0.0;
        for (std::size_t j = 0; j < slack.size(); ++j) tsl += slack[j];
        if (tsl <= zt) break;
        for (std::size_t j = 0; j < r.size(); ++j) {
            r[j] = std::max(r[j] + d * slack[j] / tsl, 0.0);
            if (capped) r[j] = std::min(r[j], xb[j]);
        }
    }
}

/** A share closure's value and Jacobian, and the joint-closure covariance. */
struct ShareValue {
    std::vector<double> s;
    Matrix<double> ds;
    std::vector<double> cn;   // Cov(S_j, N) at the means, summing to zero
    Matrix<double> dcn;       // d cn_j / d x_m, with C held fixed
};

/** How much of the second-order correction the series admits, and d tau / d ratio. */
struct ExpansionWeight {
    double tau = 1.0;
    double dtau = 0.0;
};

/**
 * Port of `local_expansion_weight` in `fluid_share_closure.m`.
 *
 * Every second-order term of the share closure is a term of the series for
 * E[1/v], whose successive terms are in the ratio `Var(v)/v^2`, so the truncation
 * is meaningful below 1 and the terms GROW above it. Nothing in the algebra
 * notices: at a near-empty station the corrections come back larger than the
 * quantity they correct, and the drift that follows is not integrable.
 *
 * One on [0,1], zero from 4 up, and the C^1 smoothstep between. Both ends matter.
 * The lower one has to be EXACTLY one on the whole convergent region, so every
 * model already inside it is bit-identical; the upper one has to be reached with a
 * vanishing derivative, because the drift is integrated and a kink in it is what
 * collapses the step size. The two thresholds are the series, not a tuning: at
 * ratio 1 successive terms stop shrinking, and at ratio 4 the standard deviation
 * of v is twice its mean, where a non-negative v has essentially no mass near the
 * point being expanded about.
 */
inline ExpansionWeight fluid_expansion_weight(double ratio) {
    const double lo = 1.0, hi = 4.0;
    ExpansionWeight w;
    if (ratio <= lo) {
        w.tau = 1.0;
        w.dtau = 0.0;
    } else if (ratio >= hi) {
        w.tau = 0.0;
        w.dtau = 0.0;
    } else {
        const double t = (ratio - lo) / (hi - lo);
        w.tau = 1.0 - t * t * (3.0 - 2.0 * t);
        w.dtau = -6.0 * t * (1.0 - t) / (hi - lo);
    }
    return w;
}

/**
 * Port of `fluid_share_closure.m`: E[w_j X_j / sum_m w_m X_m] by the delta
 * method, and its Jacobian at fixed covariance.
 *
 * THE INVARIANT TO CHECK ON ANY CHANGE HERE is that the shares sum to one. The
 * correction is exactly capacity conserving because sum_j Cov(u_j,v) = Var(v), so
 * the two correction terms cancel in the sum; a work-conserving discipline that
 * lost that identity would leak or invent capacity.
 *
 * The expansion is local and fails where v is small against its own standard
 * deviation, the exact expectation there being a Cauchy-like integral with no
 * finite mean. A raw share can come out negative; it is clipped at zero and the
 * survivors renormalised, which preserves the identity.
 */
inline ShareValue fluid_share_closure(const std::vector<double>& x, const std::vector<double>& wv,
                                      const Matrix<double>& C, bool want_jac,
                                      bool want_cov = false) {
    const std::size_t n = x.size();
    ShareValue out;
    out.s.assign(n, 0.0);
    if (want_jac) out.ds = Matrix<double>(n, n, 0.0);
    if (want_cov) {
        out.cn.assign(n, 0.0);
        if (want_jac) out.dcn = Matrix<double>(n, n, 0.0);
    }

    std::vector<double> u(n, 0.0);
    double v = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        u[j] = wv[j] * x[j];
        v += u[j];
    }
    if (v <= 0.0) return out;

    for (std::size_t j = 0; j < n; ++j) out.s[j] = u[j] / v;
    const auto plugin_jac = [&](Matrix<double>& J) {
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t m = 0; m < n; ++m)
                J(j, m) = (j == m ? wv[j] / v : 0.0) - (u[j] / (v * v)) * wv[m];
    };
    if (want_jac) plugin_jac(out.ds);

    bool have_cov = C.rows() == n && C.cols() == n;
    if (have_cov) {
        bool any = false;
        for (std::size_t a = 0; a < n && !any; ++a)
            for (std::size_t b = 0; b < n && !any; ++b)
                if (C(a, b) != 0.0) any = true;
        have_cov = any;
    }
    if (!have_cov) return out;

    std::vector<double> cuv(n, 0.0);  // Cov(u_j, v)
    double cvv = 0.0;                 // Var(v)
    for (std::size_t j = 0; j < n; ++j) {
        double acc = 0.0;
        for (std::size_t b = 0; b < n; ++b) acc += C(j, b) * wv[b];
        cuv[j] = wv[j] * acc;
        cvv += wv[j] * acc;
    }
    // How far into the series this point sits, and how much of the second-order
    // correction that leaves admissible. tau depends on x through v alone, since C
    // is held fixed here exactly as the Jacobians are.
    const ExpansionWeight tw = fluid_expansion_weight(cvv / (v * v));
    if (tw.tau <= 0.0 && !want_jac) return out;  // first-order share, already exact
    const double tau = tw.tau;
    const double dtau_scale = tw.dtau * (-2.0 * cvv / (v * v * v));  // d tau / d x_m = this * wv[m]

    if (want_cov) {
        // Cov(S_j, N) at the means, from grad(S_j)'*C*1: S_j = w_j X_j / V, so
        // dS_j/dX_a = w_j*delta_aj/V - u_j*w_a/V^2 and the two pieces contract
        // against C*1 = Cov(X,N) and w'*C*1 = Cov(V,N).
        std::vector<double> an(n, 0.0);
        double cvn = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            double acc = 0.0;
            for (std::size_t b = 0; b < n; ++b) acc += C(j, b);
            an[j] = acc;
            cvn += wv[j] * acc;
        }
        std::vector<double> cn0(n, 0.0);
        for (std::size_t j = 0; j < n; ++j) cn0[j] = wv[j] * an[j] / v - u[j] * (cvn / (v * v));
        for (std::size_t j = 0; j < n; ++j) out.cn[j] = tau * cn0[j];
        if (want_jac)
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t m = 0; m < n; ++m)
                    out.dcn(j, m) =
                        tau * (-(wv[j] * an[j]) * (wv[m] / (v * v)) -
                               (j == m ? (cvn / (v * v)) * wv[j] : 0.0) +
                               (2.0 * cvn / (v * v * v)) * u[j] * wv[m]) +
                        cn0[j] * dtau_scale * wv[m];
    }

    std::vector<double> scorr(n, 0.0);
    for (std::size_t j = 0; j < n; ++j)
        scorr[j] = -cuv[j] / (v * v) + (u[j] * cvv) / (v * v * v);
    for (std::size_t j = 0; j < n; ++j) out.s[j] += tau * scorr[j];
    if (want_jac)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t m = 0; m < n; ++m)
                out.ds(j, m) += tau * ((2.0 / (v * v * v)) * cuv[j] * wv[m] +
                                       (j == m ? (cvv / (v * v * v)) * wv[j] : 0.0) -
                                       (3.0 * cvv / (v * v * v * v)) * u[j] * wv[m]) +
                                scorr[j] * dtau_scale * wv[m];

    // A COORDINATE CARRYING NO MASS MUST NOT DECIDE THE CLIP. With u_j = 0 the
    // plug-in share is zero and the correction leaves s_j = -Cov(u_j,v)/v^2, a
    // quantity of the order of rounding whose SIGN is not meaningful. Letting it
    // select the branch below zeroes that coordinate's whole Jacobian row, and a
    // zero row is an exact zero eigenvalue: the Lyapunov step then reads the fixed
    // point as non-hyperbolic and SolverFLD silently drops from 'minnormal' to the
    // first-order method. Clip only a share that is negative BEYOND the numerical
    // zero.
    bool all_nonneg = true;
    for (std::size_t j = 0; j < n; ++j)
        all_nonneg = all_nonneg && (out.s[j] >= -lang::GlobalConstants::Zero);
    if (all_nonneg) {
        for (std::size_t j = 0; j < n; ++j)
            if (out.s[j] < 0.0) out.s[j] = 0.0;
        return out;
    }

    // Outside the region where the expansion is valid for at least one
    // coordinate: clip and renormalise the survivors so the shares still sum to
    // one. With nothing left, fall back to the plug-in share.
    std::vector<bool> act(n, false);
    double tot = 0.0;
    bool any_act = false;
    for (std::size_t j = 0; j < n; ++j) {
        act[j] = out.s[j] > 0.0;
        if (act[j]) {
            tot += out.s[j];
            any_act = true;
        }
    }
    if (!any_act) {
        for (std::size_t j = 0; j < n; ++j) out.s[j] = u[j] / v;
        if (want_jac) plugin_jac(out.ds);
        return out;
    }
    std::vector<double> snew(n, 0.0);
    if (want_jac) {
        std::vector<double> dT(n, 0.0);
        for (std::size_t m = 0; m < n; ++m)
            for (std::size_t j = 0; j < n; ++j)
                if (act[j]) dT[m] += out.ds(j, m);
        Matrix<double> dsnew(n, n, 0.0);
        for (std::size_t j = 0; j < n; ++j) {
            if (!act[j]) continue;
            for (std::size_t m = 0; m < n; ++m)
                dsnew(j, m) = out.ds(j, m) / tot - (out.s[j] / (tot * tot)) * dT[m];
        }
        out.ds = dsnew;
    }
    for (std::size_t j = 0; j < n; ++j)
        if (act[j]) snew[j] = out.s[j] / tot;
    out.s = snew;
    return out;
}

/**
 * Port of `fluid_gps_share.m`: the expected capacity share of a GPS station under
 * a normal marginal, and its Jacobian.
 *
 * GPS divides the server by WEIGHT among the BACKLOGGED classes and then equally
 * among that class's own jobs, so the share is a function of the backlog
 * INDICATOR and not of the populations. The closure is an exact ENUMERATION of
 * the 2^K patterns rather than an expansion: the share is piecewise constant over
 * them, so once the pattern probabilities are given there is no truncation error
 * left. Those probabilities come from the marginals with the continuity
 * correction P(N_k > 1/2), multiplied as if the backlogs were independent -- the
 * one approximation here, and not an innocuous one, since a closed network
 * correlates the station coordinates negatively through population conservation.
 *
 * The shares sum to 1 - P(station empty), not to 1. That is how the idle server
 * is represented: GPS is single-server, the indicator plays the role min(n,c)
 * plays at a PS station, and the caller applies no separate capacity term.
 */
inline ShareValue fluid_gps_share(const std::vector<double>& xk, const std::vector<double>& wk_in,
                                  const std::vector<double>& vk, bool want_jac) {
    const std::size_t K = xk.size();
    ShareValue out;
    out.s.assign(K, 0.0);
    if (want_jac) out.ds = Matrix<double>(K, K, 0.0);

    // The enumeration is 2^K in the classes AT ONE STATION, which is small in
    // every practical model; refuse rather than crawl.
    if (K > 12)
        throw UnsupportedError(
            "fluid_gps_share: GPS closes its capacity share by enumerating the 2^K backlog "
            "patterns of a station, and this station carries " +
            std::to_string(K) +
            " classes. Above 12 the enumeration is no longer tractable; use method 'closing' with a "
            "DPS station instead");

    double sw = 0.0;
    for (std::size_t k = 0; k < K; ++k) sw += wk_in[k];
    if (sw <= 0.0) return out;
    std::vector<double> wk(K, 0.0);
    for (std::size_t k = 0; k < K; ++k) wk[k] = wk_in[k] / sw;

    // P(class k backlogged) = P(N_k >= 1) for an INTEGER population, so the
    // normal approximation takes the continuity correction P(N_k > 1/2).
    // Thresholding at 1 instead understates the backlog probability, and at
    // sigma = 0 it makes the share vanish for every class with x_k < 1, which
    // stalls the server completely and is an ABSORBING state for the ODE. The
    // sigma = 0 fallback is the fluid statement that positive mass is backlogged.
    std::vector<double> p(K, 0.0), dp(K, 0.0);
    for (std::size_t k = 0; k < K; ++k) {
        if (vk[k] > 0.0) {
            const double sd = std::sqrt(vk[k]);
            const double z = (xk[k] - 0.5) / sd;
            p[k] = closure_normcdf(z);
            dp[k] = closure_normpdf(z) / sd;
        } else {
            p[k] = (xk[k] > 0.0) ? 1.0 : 0.0;
            dp[k] = 0.0;
        }
    }

    Matrix<double> dsdp(K, K, 0.0);
    const unsigned long long masks = 1ULL << K;
    for (unsigned long long mask = 1; mask < masks; ++mask) {
        double W = 0.0;
        for (std::size_t k = 0; k < K; ++k)
            if (mask & (1ULL << k)) W += wk[k];
        if (W <= 0.0) continue;  // every backlogged class here carries zero weight
        std::vector<double> q(K, 0.0);
        double prodq = 1.0;
        for (std::size_t k = 0; k < K; ++k) {
            q[k] = (mask & (1ULL << k)) ? p[k] : (1.0 - p[k]);
            prodq *= q[k];
        }
        for (std::size_t k = 0; k < K; ++k)
            if (mask & (1ULL << k)) out.s[k] += prodq * wk[k] / W;
        if (!want_jac) continue;
        for (std::size_t m = 0; m < K; ++m) {
            double prodm = 1.0;  // the product over j != m
            for (std::size_t j = 0; j < K; ++j)
                if (j != m) prodm *= q[j];
            const double sgn = (mask & (1ULL << m)) ? 1.0 : -1.0;
            for (std::size_t k = 0; k < K; ++k)
                if (mask & (1ULL << k)) dsdp(k, m) += sgn * prodm * wk[k] / W;
        }
    }

    if (want_jac)
        for (std::size_t k = 0; k < K; ++k)
            for (std::size_t m = 0; m < K; ++m) out.ds(k, m) = dsdp(k, m) * dp[m];
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_CLOSURES_H
