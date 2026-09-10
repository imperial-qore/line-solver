/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SJN_H
#define LINE_API_PFQN_SJN_H

/**
 * Closed networks with non-preemptive shortest-job-next (SJN/SJF) stations.
 *
 * Port of matlab/src/api/pfqn/pfqn_mvasjn.m and pfqn_amvasjn.m together with
 * the six private helpers they share (sjn_args, sjn_fit, sjn_setup, sjn_quad,
 * sjn_station, sjn_cap), cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/{Pfqn_mvasjn,Pfqn_amvasjn,SjnSupport}.java.
 *
 * THE EQUATION. A tagged customer whose service requirement is x waits, by the
 * arrival theorem, for the residual life of the job in service, for the work of
 * the queued jobs that will be served before it, and for the work of the jobs
 * that overtake it while it waits (Kant 1992, eqs. 1-7):
 *
 *   W(x,n) = [ (1+CV^2) s U(n-1)/2 + X(n-1) phi(x,n-1) ] / [ 1 - X(n-1) theta(x) ]
 *   theta(x) = int_0^x t f(t) dt,   phi(x,n) = int_0^x W(t,n) t f(t) dt
 *   R(n) = s + int_0^inf W(x,n) f(x) dx
 *
 * W(.,n) needs only phi(.,n-1), so the profile is carried alongside the
 * population recursion. `pfqn_mvasjn` does exactly that over the whole lattice;
 * `pfqn_amvasjn` replaces the lattice by a Schweitzer closure applied to the
 * SIZE-RESOLVED queue length lam_k W_k(x) f_k(x) rather than to its integral,
 * which is why the two share `sjn_station` unchanged and differ only in the
 * deflation vector beta handed to it.
 *
 * THE SIZE DENSITY IS NOT AN INPUT. Only its mean and SCV are, and the density
 * is reconstructed by the two-moment branching-Erlang fit the reference
 * prescribes. That is what makes theta and the tail integrals closed form and
 * lets the x-integrals run on a FIXED grid: W(.,n) is needed again at the next
 * population step, so a rule that samples at arbitrary abscissae cannot be
 * used. Beyond the grid edge Lx the profile is closed by the analytic tail
 * W(x,n) = a - b exp(-c (x - Lx)) of eqs. (11)-(14).
 *
 * ARITHMETIC: DOUBLE, NOT TEMPLATED. Unlike its pfqn neighbours this header is
 * not generic in the number type. The recursion evaluates the REGULARIZED
 * INCOMPLETE GAMMA in both its branches, for which the port has no T-generic
 * implementation, and its accuracy is set by a 33-point Simpson grid rather
 * than by the arithmetic, so a wider T would buy nothing. Callers gate on
 * num_traits<T>::has_transcendental and convert at the boundary; the exact
 * (Rational) path refuses SJN by name in solver_mva_sjn.h.
 *
 * WHY NOT REUSE mam::gammainc_lower. It returns only P(a,x). The tail
 * correction needs Q(a,x) at magnitudes around 1e-300, where 1 - P is exactly
 * zero, so an upper branch that is computed rather than subtracted is
 * mandatory. `detail::gammainc` below returns both from the one continued
 * fraction, which is also what MATLAB's gammainc(...,'upper') does.
 *
 * WARNINGS BECOME FLAGS. The reference calls line_warning when the utilization
 * cap binds and when the fixed point runs out of iterations. This port has no
 * warning channel, so `SjnResult::capped` and `SjnResult::converged` carry the
 * same information to the caller. The hard `line_error` cases stay exceptions.
 *
 * Reference: K. Kant, "MVA approximations for SJN scheduling", Performance
 * Evaluation 15(1):41-61, 1992.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Options of the SJN solvers, the fields sjn_args fills in. */
struct SjnOptions {
    std::size_t ns = 32;        ///< grid subdivisions, must be even
    double Lfactor = 8.0;       ///< grid extent in units of the largest mean service time
    std::vector<int> prio;      ///< distinct levels, 1 = highest; empty for the pooled reading
    double tol = 1e-8;
    std::size_t iter_max = 1000;
    double umax = 0.999;        ///< utilization cap, strictly below one
};

/** The conditional waiting time profile at one SJN station, the reference's WX. */
struct SjnProfile {
    std::size_t station = 0;   ///< 1-based row index into L
    std::vector<double> x;     ///< the grid
    Matrix<double> W;          ///< (ngrid x R)
    Matrix<double> tail;       ///< (R x 3), the tail parameters a, b, c
};

/** Return block of pfqn_mvasjn and pfqn_amvasjn. */
struct SjnResult {
    std::vector<double> XN;          ///< (R)
    Matrix<double> QN, UN, CN;       ///< (M x R)
    std::vector<SjnProfile> WX;
    std::size_t iter = 1;
    bool converged = true;           ///< always true for the lattice recursion
    bool capped = false;             ///< the utilization cap was binding somewhere
};

/**
 * The conditional waiting time equation has no solution at some population.
 *
 * Raised by pfqn_mvasjn only: the cap rescales the profile that the next
 * population step reads back, so the correction compounds along the lattice
 * and the recursion oscillates. pfqn_amvasjn solves the same profile by a
 * fixed point and does cap, the iteration being self-consistent.
 */
class SjnStarvationError : public NumericError {
public:
    explicit SjnStarvationError(const std::string& what) : NumericError(what) {}
};

namespace detail {

/**
 * Regularized incomplete gamma, both branches at once.
 *
 * Series below the transition point, continued fraction above it (Numerical
 * Recipes 6.2). The branch that the chosen expansion computes DIRECTLY is
 * returned without a subtraction, so the small one keeps its relative accuracy
 * down to the underflow limit, which the tail integrals depend on.
 */
struct GammaInc {
    double p;  ///< lower P(a,x)
    double q;  ///< upper Q(a,x)
};

inline GammaInc gammainc(double a, double x) {
    if (a <= 0.0) throw InputError("pfqn_sjn: the incomplete gamma needs a positive order");
    if (x < 0.0) throw InputError("pfqn_sjn: the incomplete gamma needs a non-negative argument");
    if (x == 0.0) return GammaInc{0.0, 1.0};
    const double gln = std::lgamma(a);
    const double lead = std::exp(-x + a * std::log(x) - gln);
    if (x < a + 1.0) {
        double ap = a, del = 1.0 / a, sum = del;
        for (int n = 0; n < 1000; ++n) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
        }
        const double p = sum * lead;
        return GammaInc{p, 1.0 - p};
    }
    const double tiny = 1e-300;
    double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
    for (int i = 1; i <= 1000; ++i) {
        const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::fabs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-16) break;
    }
    const double q = lead * h;
    return GammaInc{1.0 - q, q};
}

/** Erlang mixture fitted to a mean and an SCV, the reference's sjn_fit. */
struct SjnFit {
    std::vector<double> w;
    std::vector<int> k;
    std::vector<double> mu;
    bool empty() const { return w.empty(); }
};

/**
 * Branching Erlang (Erlang(k-1) and Erlang(k) sharing a rate) below CV^2 = 1,
 * balanced-means hyperexponential above it, exponential at it.
 */
inline SjnFit sjn_fit(double s, double cv2) {
    SjnFit f;
    if (s <= 0.0) return f;
    if (cv2 < 0.0) throw InputError("pfqn_sjn: negative squared coefficient of variation");
    if (std::fabs(cv2 - 1.0) < 1e-8) {
        f.w = {1.0};
        f.k = {1};
        f.mu = {1.0 / s};
    } else if (cv2 < 1.0) {
        const double kd = std::ceil(1.0 / cv2);
        const double p = (kd * cv2 - std::sqrt(kd * (1.0 + cv2) - kd * kd * cv2)) / (1.0 + cv2);
        const double mu = (kd - p) / s;
        f.w = {p, 1.0 - p};
        f.k = {static_cast<int>(kd) - 1, static_cast<int>(kd)};
        f.mu = {mu, mu};
    } else {
        const double p = 0.5 * (1.0 + std::sqrt((cv2 - 1.0) / (cv2 + 1.0)));
        f.w = {p, 1.0 - p};
        f.k = {1, 1};
        f.mu = {2.0 * p / s, 2.0 * (1.0 - p) / s};
    }
    return f;
}

/** Density of the mixture on a grid, the reference's sjn_quad('pdf'). */
inline std::vector<double> sjn_pdf(const SjnFit& f, const std::vector<double>& x) {
    // realmin, not the smallest denormal: MATLAB's max(x,realmin) is the model.
    const double floorx = std::numeric_limits<double>::min();
    std::vector<double> y(x.size(), 0.0);
    for (std::size_t j = 0; j < f.w.size(); ++j) {
        const double kd = static_cast<double>(f.k[j]), mu = f.mu[j];
        for (std::size_t i = 0; i < x.size(); ++i)
            y[i] += f.w[j] * std::exp(kd * std::log(mu) +
                                      (kd - 1.0) * std::log(std::max(x[i], floorx)) - mu * x[i] -
                                      std::lgamma(kd));
    }
    return y;
}

/** int_0^x t f(t) dt in closed form, the reference's sjn_quad('theta'). */
inline std::vector<double> sjn_theta(const SjnFit& f, const std::vector<double>& x) {
    std::vector<double> y(x.size(), 0.0);
    for (std::size_t j = 0; j < f.w.size(); ++j) {
        const double kd = static_cast<double>(f.k[j]), mu = f.mu[j];
        for (std::size_t i = 0; i < x.size(); ++i)
            y[i] += f.w[j] * (kd / mu) * gammainc(kd + 1.0, mu * x[i]).p;
    }
    return y;
}

/** int_x^inf f(t) dt in closed form, the reference's sjn_quad('ccdf'). */
inline double sjn_ccdf(const SjnFit& f, double x) {
    double y = 0.0;
    for (std::size_t j = 0; j < f.w.size(); ++j)
        y += f.w[j] * gammainc(static_cast<double>(f.k[j]), f.mu[j] * x).q;
    return y;
}

/**
 * int_Lx^inf t^order exp(-c (t-Lx)) f(t) dt, the reference's sjn_quad('tailmom').
 * Evaluated in logarithms so that exp(c Lx) cannot overflow against an
 * underflowing incomplete gamma.
 */
inline double sjn_tailmom(const SjnFit& f, double Lx, double c, int order) {
    double y = 0.0;
    for (std::size_t j = 0; j < f.w.size(); ++j) {
        const double kd = static_cast<double>(f.k[j]), mu = f.mu[j];
        const double rate = mu + c;
        const double g = gammainc(kd + static_cast<double>(order), rate * Lx).q;
        if (g <= 0.0) continue;
        double lg = c * Lx + kd * std::log(mu / rate) + std::log(g);
        if (order == 1) lg += std::log(kd / rate);
        y += f.w[j] * std::exp(lg);
    }
    return y;
}

/** Composite Simpson over an even number of subdivisions. */
inline double sjn_simpson(const std::vector<double>& y, double dx) {
    const std::size_t n = y.size();
    double odd = 0.0, even = 0.0;
    for (std::size_t i = 1; i + 1 < n; i += 2) odd += y[i];
    for (std::size_t i = 2; i + 1 < n; i += 2) even += y[i];
    return dx / 3.0 * (y[0] + y[n - 1] + 4.0 * odd + 2.0 * even);
}

/**
 * Cumulative Simpson: full panels at the even nodes and a half panel at the odd
 * ones, so the primitive is available at EVERY grid node. Quadrature at
 * arbitrary abscissae could not provide that, the profile being needed again at
 * the next population step.
 */
inline std::vector<double> sjn_cumsimpson(const std::vector<double>& y, double dx) {
    const std::size_t n = y.size();
    std::vector<double> I(n, 0.0);
    for (std::size_t i = 2; i < n; i += 2)
        I[i] = I[i - 2] + dx / 3.0 * (y[i - 2] + 4.0 * y[i - 1] + y[i]);
    for (std::size_t i = 1; i < n; i += 2) {
        if (i + 1 < n)
            I[i] = I[i - 1] + dx / 12.0 * (5.0 * y[i - 1] + 8.0 * y[i] - y[i + 1]);
        else
            I[i] = I[i - 1] + dx / 12.0 * (-y[i - 2] + 8.0 * y[i - 1] + 5.0 * y[i]);
    }
    return I;
}

/** Grid of one SJN station and the population-independent integrals over it. */
struct SjnGrid {
    double Lx = 0.0, dx = 0.0;
    std::vector<double> x;
    Matrix<double> f, theta;   ///< (ngrid x R)
    std::vector<double> tail0, tail1;
    std::vector<SjnFit> fit;
};

/** Port of sjn_setup. */
inline SjnGrid sjn_setup(const std::vector<double>& S, const std::vector<double>& scv,
                         std::size_t ns, double Lfactor) {
    const std::size_t R = S.size(), ngrid = ns + 1;
    double smax = 0.0;
    for (double s : S) smax = std::max(smax, s);
    if (smax <= 0.0)
        throw InputError("pfqn_sjn: the station has zero service demand in every class");
    SjnGrid G;
    G.Lx = Lfactor * smax;
    G.dx = G.Lx / static_cast<double>(ns);
    G.x.assign(ngrid, 0.0);
    for (std::size_t i = 0; i < ngrid; ++i) G.x[i] = static_cast<double>(i) * G.dx;
    G.f = Matrix<double>(ngrid, R, 0.0);
    G.theta = Matrix<double>(ngrid, R, 0.0);
    G.tail0.assign(R, 0.0);
    G.tail1.assign(R, 0.0);
    G.fit.assign(R, SjnFit());
    for (std::size_t r = 0; r < R; ++r) {
        G.fit[r] = sjn_fit(S[r], scv[r]);
        if (G.fit[r].empty()) continue;
        const std::vector<double> fr = sjn_pdf(G.fit[r], G.x);
        const std::vector<double> tr = sjn_theta(G.fit[r], G.x);
        for (std::size_t i = 0; i < ngrid; ++i) {
            G.f(i, r) = fr[i];
            G.theta(i, r) = tr[i];
        }
        G.tail0[r] = sjn_ccdf(G.fit[r], G.Lx);
        G.tail1[r] = S[r] - tr[ngrid - 1];
    }
    return G;
}

/** State of one SJN station at the reference population, the reference's `st`. */
struct SjnState {
    std::vector<double> lam, U, Q;
    const Matrix<double>* W = nullptr;    ///< (ngrid x R)
    const Matrix<double>* phi = nullptr;  ///< (ngrid x R)
    std::vector<double> phiinf;
};

/** Outcome of one evaluation of the conditional waiting time equation. */
struct SjnStationResult {
    double C = 0.0;
    std::vector<double> W, phi;
    double phiinf = 0.0;
    double tail[3] = {0.0, 0.0, 0.0};
};

/** `m` is 0-based here; the message reports the reference's 1-based index. */
inline NumericError sjn_singular(std::size_t m) {
    return NumericError(
        "pfqn_sjn: the SJN recursion at station " + std::to_string(m + 1) +
        " has no solution: the work brought by jobs no longer than the tagged one saturates the "
        "server, at which point long jobs starve and the arrival theorem no longer holds. Reduce "
        "the load at that station or model it with SolverCTMC or SolverLDES");
}

/**
 * Port of sjn_station: one evaluation of the SJN conditional waiting time
 * equation for a tagged customer of class r at station m.
 *
 * lam_k W_k(x) f_k(x) is the density, in the job size x, of the queued class-k
 * customers, so deflating it by beta_k turns the same equation into either the
 * exact recursion (beta = 1, the state already being the one at n - e_r) or the
 * Schweitzer closure (beta_r = (N_r-1)/N_r, the state being the one at N).
 */
inline SjnStationResult sjn_station(std::size_t m, std::size_t r, const SjnGrid& G,
                                    const std::vector<double>& S, const std::vector<double>& scv,
                                    const std::vector<double>& V, const SjnState& st,
                                    const std::vector<double>& beta, bool useprio,
                                    const std::vector<int>& prio) {
    const std::size_t R = S.size(), ngrid = G.x.size();
    std::vector<double> lamb(R), Ub(R), Qb(R);
    double RL = 0.0;
    for (std::size_t k = 0; k < R; ++k) {
        lamb[k] = beta[k] * st.lam[k];
        Ub[k] = beta[k] * st.U[k];
        Qb[k] = beta[k] * st.Q[k];
        RL += (1.0 + scv[k]) * S[k] * Ub[k] / 2.0;
    }
    std::vector<double> num(ngrid), den(ngrid);
    double numinf = 0.0, deninf = 0.0;
    if (useprio) {
        double base = RL, uhi = 0.0;
        for (std::size_t k = 0; k < R; ++k)
            if (prio[k] < prio[r]) {
                base += S[k] * (Qb[k] - Ub[k]);
                uhi += Ub[k];
            }
        for (std::size_t i = 0; i < ngrid; ++i) {
            num[i] = base + lamb[r] * (*st.phi)(i, r);
            den[i] = 1.0 - uhi - lamb[r] * G.theta(i, r);
        }
        numinf = base + lamb[r] * st.phiinf[r];
        deninf = 1.0 - uhi - lamb[r] * S[r];
    } else {
        double usum = 0.0, phiinfsum = 0.0;
        for (std::size_t k = 0; k < R; ++k) {
            usum += lamb[k] * S[k];
            phiinfsum += lamb[k] * st.phiinf[k];
        }
        for (std::size_t i = 0; i < ngrid; ++i) {
            double n = RL, d = 1.0;
            for (std::size_t k = 0; k < R; ++k) {
                n += lamb[k] * (*st.phi)(i, k);
                d -= lamb[k] * G.theta(i, k);
            }
            num[i] = n;
            den[i] = d;
        }
        numinf = RL + phiinfsum;
        deninf = 1.0 - usum;
    }
    for (std::size_t i = 0; i < ngrid; ++i)
        if (den[i] <= 0.0) throw sjn_singular(m);
    if (deninf <= 0.0) throw sjn_singular(m);

    SjnStationResult out;
    out.W.assign(ngrid, 0.0);
    for (std::size_t i = 0; i < ngrid; ++i) out.W[i] = num[i] / den[i];
    const double Winf = numinf / deninf;
    // eq. (14), generalised by differentiating the recursion at the grid edge
    double slope;
    if (useprio) {
        slope = G.Lx * lamb[r] * G.f(ngrid - 1, r) *
                ((*st.W)(ngrid - 1, r) + out.W[ngrid - 1]) / den[ngrid - 1];
    } else {
        double acc = 0.0;
        for (std::size_t k = 0; k < R; ++k)
            acc += lamb[k] * G.f(ngrid - 1, k) * ((*st.W)(ngrid - 1, k) + out.W[ngrid - 1]);
        slope = G.Lx * acc / den[ngrid - 1];
    }
    double a = Winf, b = Winf - out.W[ngrid - 1], c;
    if (b <= 0.0) {
        b = 0.0;
        c = 0.0;
    } else if (slope < 0.0) {
        throw NumericError("pfqn_sjn: the conditional waiting time at SJN station " +
                           std::to_string(m + 1) +
                           " decreases in the job size, which the discipline forbids: the "
                           "recursion has become numerically unstable");
    } else {
        c = slope / b;
    }
    out.tail[0] = a;
    out.tail[1] = b;
    out.tail[2] = c;

    std::vector<double> integrand(ngrid);
    for (std::size_t i = 0; i < ngrid; ++i) integrand[i] = out.W[i] * G.x[i] * G.f(i, r);
    out.phi = sjn_cumsimpson(integrand, G.dx);
    out.phiinf = out.phi[ngrid - 1] + a * G.tail1[r] - b * sjn_tailmom(G.fit[r], G.Lx, c, 1);
    for (std::size_t i = 0; i < ngrid; ++i) integrand[i] = out.W[i] * G.f(i, r);
    const double Wbar =
        sjn_simpson(integrand, G.dx) + a * G.tail0[r] - b * sjn_tailmom(G.fit[r], G.Lx, c, 0);
    out.C = V[r] * (S[r] + Wbar);
    return out;
}

/** Throughputs implied by the residence times, keeping Little's law exact. */
inline std::vector<double> sjn_thru(const Matrix<double>& C, const std::vector<double>& N,
                                    const std::vector<double>& Z) {
    const std::size_t M = C.rows(), R = N.size();
    std::vector<double> X(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] <= 0.0) continue;
        double den = Z[r];
        for (std::size_t m = 0; m < M; ++m) den += C(m, r);
        X[r] = N[r] / den;
    }
    return X;
}

/** Outcome of the utilization cap. */
struct SjnCapResult {
    std::vector<double> X, kappa;
    bool bound = false;
};

/**
 * Port of sjn_cap: enforce U <= umax at every SJN station by inflating its
 * waiting time, and return the throughputs the capped residence times imply.
 *
 * The SJN response time equation is an open-system one: its denominator is
 * 1 - U(x), and it has no solution once that reaches one. A closed network
 * never reaches it in reality, but the approximation can, because it
 * underestimates the residence time at a congested SJN station and the
 * resulting throughput then exceeds the station capacity. What is imposed is
 * the utilization law sum_r X_r L_mr <= umax, an exact property of the network
 * and not a property of the approximation. It acts on the EXCESS C - L, never
 * on the throughput, so that X (Z + sum_m C) = N still holds exactly and no
 * jobs are lost.
 */
inline SjnCapResult sjn_cap(Matrix<double>& C, const Matrix<double>& L,
                            const std::vector<double>& N, const std::vector<double>& Z,
                            const std::vector<std::size_t>& sjnset, double umax) {
    const std::size_t nsjn = sjnset.size(), R = N.size();
    SjnCapResult out;
    out.kappa.assign(nsjn, 1.0);
    out.X = sjn_thru(C, N, Z);
    if (nsjn == 0) return out;
    if (umax >= 1.0)
        throw InputError("pfqn_sjn: the utilization cap must be strictly below one, the response "
                         "time equation is singular at one");
    // rho at station m had the excess been inflated by kappa, C left untouched
    const auto rho_at = [&](std::size_t m, const std::vector<double>& Wq, double kappa) {
        std::vector<double> saved(R);
        for (std::size_t r = 0; r < R; ++r) {
            saved[r] = C(m, r);
            C(m, r) = L(m, r) + kappa * Wq[r];
        }
        const std::vector<double> X = sjn_thru(C, N, Z);
        double rho = 0.0;
        for (std::size_t r = 0; r < R; ++r) rho += X[r] * L(m, r);
        for (std::size_t r = 0; r < R; ++r) C(m, r) = saved[r];
        return rho;
    };
    for (int sweep = 0; sweep < 20; ++sweep) {
        bool viol = false;
        for (std::size_t q = 0; q < nsjn; ++q) {
            const std::size_t m = sjnset[q];
            double rho = 0.0;
            for (std::size_t r = 0; r < R; ++r) rho += out.X[r] * L(m, r);
            if (rho <= umax) continue;
            viol = true;
            out.bound = true;
            std::vector<double> Wq(R);
            for (std::size_t r = 0; r < R; ++r) Wq[r] = C(m, r) - L(m, r);
            double hi = 2.0;
            while (rho_at(m, Wq, hi) > umax) {
                hi *= 2.0;
                if (hi > 1e12)
                    throw NumericError(
                        "pfqn_sjn: station " + std::to_string(m + 1) +
                        " cannot be brought under the utilization cap by any waiting time: its "
                        "service demands alone saturate it at this population");
            }
            double lo = 1.0;
            for (int b = 0; b < 200; ++b) {
                const double mid = (lo + hi) / 2.0;
                if (rho_at(m, Wq, mid) > umax)
                    lo = mid;
                else
                    hi = mid;
            }
            out.kappa[q] *= hi;
            for (std::size_t r = 0; r < R; ++r) C(m, r) = L(m, r) + hi * Wq[r];
            out.X = sjn_thru(C, N, Z);
        }
        if (!viol) return out;
    }
    throw NumericError("pfqn_sjn: the utilization cap did not settle across the SJN stations");
}

/** Normalised arguments, the reference's sjn_args. */
struct SjnArgs {
    std::size_t M = 0, R = 0;
    std::vector<double> N, Z;
    Matrix<double> scv, V, S;
    std::vector<std::size_t> sjnset;  ///< 0-based rows of L
    SjnOptions opt;
};

/** Port of sjn_args. */
inline SjnArgs sjn_args(const Matrix<double>& L, const std::vector<double>& N,
                        const std::vector<double>& Z, const Matrix<double>& scv,
                        const std::vector<std::size_t>& sjnset, const Matrix<double>& V,
                        const SjnOptions& options) {
    SjnArgs a;
    a.M = L.rows();
    a.R = L.cols();
    a.opt = options;
    if (N.size() != a.R)
        throw InputError("pfqn_sjn: demand matrix and population vector have different number of "
                         "classes");
    a.N.assign(a.R, 0.0);
    for (std::size_t r = 0; r < a.R; ++r) {
        a.N[r] = std::round(N[r]);
        if (a.N[r] < 0.0) throw InputError("pfqn_sjn: negative class populations");
    }
    a.Z = Z.empty() ? std::vector<double>(a.R, 0.0) : Z;
    if (a.Z.size() != a.R) throw InputError("pfqn_sjn: Z has the wrong length");
    a.scv = scv.empty() ? Matrix<double>(a.M, a.R, 1.0) : scv;
    if (a.scv.rows() != a.M || a.scv.cols() != a.R)
        throw InputError("pfqn_sjn: scv has the wrong shape");
    a.sjnset = sjnset;
    for (std::size_t i = 0; i < a.sjnset.size(); ++i) {
        if (a.sjnset[i] >= a.M)
            throw InputError("pfqn_sjn: sjnset contains a station index outside 1..M");
        for (std::size_t j = i + 1; j < a.sjnset.size(); ++j)
            if (a.sjnset[i] == a.sjnset[j])
                throw InputError("pfqn_sjn: sjnset repeats a station index");
    }
    a.V = V.empty() ? Matrix<double>(a.M, a.R, 1.0) : V;
    if (a.V.rows() != a.M || a.V.cols() != a.R)
        throw InputError("pfqn_sjn: V has the wrong shape");
    // the job size the discipline compares is one VISIT's service time, not the
    // demand accumulated over all visits
    a.S = Matrix<double>(a.M, a.R, 0.0);
    for (std::size_t m = 0; m < a.M; ++m)
        for (std::size_t r = 0; r < a.R; ++r)
            if (a.V(m, r) > 0.0) a.S(m, r) = L(m, r) / a.V(m, r);
    if (a.opt.ns == 0 || a.opt.ns % 2 != 0)
        throw InputError("pfqn_sjn: options.ns must be even, composite Simpson integrates over "
                         "panels of two subdivisions");
    if (a.opt.umax <= 0.0 || a.opt.umax >= 1.0)
        throw InputError("pfqn_sjn: options.umax must lie strictly between zero and one");
    if (!a.opt.prio.empty()) {
        if (a.opt.prio.size() != a.R)
            throw InputError("pfqn_sjn: options.prio must have one priority level per class");
        for (std::size_t i = 0; i < a.R; ++i)
            for (std::size_t j = i + 1; j < a.R; ++j)
                if (a.opt.prio[i] == a.opt.prio[j])
                    throw InputError("pfqn_sjn: options.prio must assign distinct levels, ties "
                                     "across classes are not covered by the SJN priority "
                                     "equations");
    }
    return a;
}

}  // namespace detail

/**
 * Exact-lattice MVA for closed networks with SJN stations, the unidirectional
 * scheme of Kant 1992.
 *
 * The recursion is explicit -- W(.,n) needs only phi(.,n-1) -- so the profile
 * is carried alongside the population recursion and the whole lattice of
 * prod(N+1) states is stepped through. Two multiclass readings of "shortest
 * job" are selected by options.prio: POOLED (empty, the default) compares the
 * jobs of every class by size directly and collapses to eq. (6) for a single
 * class; PRIORITY (distinct levels, 1 = highest) is method A of eq. (21), SJN
 * applying only within a class. Method B, which evaluates the denominator at
 * the non-integral population n - Q(n), is not implemented in either codebase.
 *
 * @param L      (M x R) demands at the queueing stations
 * @param N      (R) populations
 * @param Z      (R) think times, empty for none
 * @param scv    (M x R) squared coefficients of variation, empty for ones
 * @param sjnset 0-based rows of L that schedule by SJN
 * @param V      (M x R) visit ratios, empty for ones
 * @param options quadrature grid, tolerances and iteration caps of the recursion
 */
inline SjnResult pfqn_mvasjn(const Matrix<double>& L, const std::vector<double>& N,
                             const std::vector<double>& Z, const Matrix<double>& scv,
                             const std::vector<std::size_t>& sjnset, const Matrix<double>& V,
                             const SjnOptions& options) {
    const detail::SjnArgs a = detail::sjn_args(L, N, Z, scv, sjnset, V, options);
    const std::size_t M = a.M, R = a.R, nsjn = a.sjnset.size();
    const bool useprio = !a.opt.prio.empty();
    const std::size_t ns = a.opt.ns, ngrid = ns + 1;

    std::vector<detail::SjnGrid> G(nsjn);
    for (std::size_t q = 0; q < nsjn; ++q) {
        std::vector<double> Sq(R), scvq(R);
        for (std::size_t r = 0; r < R; ++r) {
            Sq[r] = a.S(a.sjnset[q], r);
            scvq[r] = a.scv(a.sjnset[q], r);
        }
        G[q] = detail::sjn_setup(Sq, scvq, ns, a.opt.Lfactor);
    }

    // MATLAB reaches the same limit as an out-of-memory on prod(N+1); say so
    // rather than wrap the index arithmetic silently
    double lattice = 1.0;
    for (std::size_t r = 0; r < R; ++r) lattice *= a.N[r] + 1.0;
    if (lattice > 1e9)
        throw InputError("pfqn_mvasjn: the population lattice has " + std::to_string(lattice) +
                         " states and does not fit; use pfqn_amvasjn (method 'amva')");
    std::vector<std::size_t> stride(R, 1);
    std::size_t npop = 1;
    for (std::size_t r = 0; r < R; ++r) {
        stride[r] = npop;
        npop *= static_cast<std::size_t>(a.N[r]) + 1;
    }

    // one slab per lattice point: X, [Q,U,C] per station, and per SJN station
    // the profile W, its primitive phi, phi at infinity and the tail parameters
    std::vector<double> Xp(npop * R, 0.0), Qp(npop * M * R, 0.0), Up(npop * M * R, 0.0),
        Cp(npop * M * R, 0.0);
    std::vector<Matrix<double>> Wp(nsjn * npop), Pp(nsjn * npop);
    std::vector<double> Ip(nsjn * npop * R, 0.0), Tp(nsjn * npop * R * 3, 0.0);
    for (auto& m : Wp) m = Matrix<double>(ngrid, R, 0.0);
    for (auto& m : Pp) m = Matrix<double>(ngrid, R, 0.0);

    SjnResult res;
    Matrix<double> Call(M, R, 0.0);
    for (std::size_t idx = 1; idx < npop; ++idx) {
        std::vector<double> n(R, 0.0);
        for (std::size_t r = 0; r < R; ++r)
            n[r] = static_cast<double>((idx / stride[r]) %
                                       (static_cast<std::size_t>(a.N[r]) + 1));
        Call.fill(0.0);
        for (std::size_t r = 0; r < R; ++r) {
            if (n[r] == 0.0) continue;
            const std::size_t iprev = idx - stride[r];
            for (std::size_t m = 0; m < M; ++m) {
                std::size_t q = nsjn;
                for (std::size_t t = 0; t < nsjn; ++t)
                    if (a.sjnset[t] == m) q = t;
                if (q == nsjn) {
                    double qsum = 0.0;
                    for (std::size_t k = 0; k < R; ++k) qsum += Qp[(iprev * M + m) * R + k];
                    Call(m, r) = L(m, r) * (1.0 + qsum);
                    continue;
                }
                // the population step already supplies the neighbouring profile,
                // so no deflation is needed
                const std::vector<double> beta(R, 1.0);
                detail::SjnState st;
                st.lam.assign(R, 0.0);
                st.U.assign(R, 0.0);
                st.Q.assign(R, 0.0);
                st.phiinf.assign(R, 0.0);
                for (std::size_t k = 0; k < R; ++k) {
                    st.lam[k] = Xp[iprev * R + k] * a.V(m, k);
                    st.U[k] = Up[(iprev * M + m) * R + k];
                    st.Q[k] = Qp[(iprev * M + m) * R + k];
                    st.phiinf[k] = Ip[(q * npop + iprev) * R + k];
                }
                st.W = &Wp[q * npop + iprev];
                st.phi = &Pp[q * npop + iprev];
                std::vector<double> Sm(R), scvm(R), Vm(R);
                for (std::size_t k = 0; k < R; ++k) {
                    Sm[k] = a.S(m, k);
                    scvm[k] = a.scv(m, k);
                    Vm[k] = a.V(m, k);
                }
                const detail::SjnStationResult sr =
                    detail::sjn_station(m, r, G[q], Sm, scvm, Vm, st, beta, useprio, a.opt.prio);
                Call(m, r) = sr.C;
                for (std::size_t i = 0; i < ngrid; ++i) {
                    Wp[q * npop + idx](i, r) = sr.W[i];
                    Pp[q * npop + idx](i, r) = sr.phi[i];
                }
                Ip[(q * npop + idx) * R + r] = sr.phiinf;
                for (int t = 0; t < 3; ++t)
                    Tp[((q * npop + idx) * R + r) * 3 + t] = sr.tail[t];
            }
        }
        const detail::SjnCapResult cap = detail::sjn_cap(Call, L, n, a.Z, a.sjnset, a.opt.umax);
        if (cap.bound) {
            // the cap has invalidated the profile the next population step reads back
            std::string npos;
            for (std::size_t r = 0; r < R; ++r)
                npos += (r ? " " : "") + std::to_string(static_cast<long long>(n[r]));
            throw SjnStarvationError(
                "pfqn_mvasjn: the utilization cap of " + std::to_string(a.opt.umax) +
                " was binding at an SJN station at population [" + npos +
                "]: the station is in the starvation regime, where the conditional waiting time "
                "equation has no solution and the population lattice no valid continuation. Use "
                "the Schweitzer fixed point (pfqn_amvasjn, method 'amva'), SolverCTMC or "
                "SolverLDES.");
        }
        for (std::size_t r = 0; r < R; ++r) Xp[idx * R + r] = cap.X[r];
        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t r = 0; r < R; ++r) {
                Cp[(idx * M + m) * R + r] = Call(m, r);
                Qp[(idx * M + m) * R + r] = cap.X[r] * Call(m, r);
                Up[(idx * M + m) * R + r] = cap.X[r] * L(m, r);
            }
    }

    const std::size_t last = npop - 1;
    res.XN.assign(R, 0.0);
    res.QN = Matrix<double>(M, R, 0.0);
    res.UN = Matrix<double>(M, R, 0.0);
    res.CN = Matrix<double>(M, R, 0.0);
    for (std::size_t r = 0; r < R; ++r) res.XN[r] = Xp[last * R + r];
    for (std::size_t m = 0; m < M; ++m)
        for (std::size_t r = 0; r < R; ++r) {
            res.QN(m, r) = Qp[(last * M + m) * R + r];
            res.UN(m, r) = Up[(last * M + m) * R + r];
            res.CN(m, r) = Cp[(last * M + m) * R + r];
        }
    res.WX.resize(nsjn);
    for (std::size_t q = 0; q < nsjn; ++q) {
        res.WX[q].station = a.sjnset[q] + 1;
        res.WX[q].x = G[q].x;
        res.WX[q].W = Wp[q * npop + last];
        res.WX[q].tail = Matrix<double>(R, 3, 0.0);
        for (std::size_t r = 0; r < R; ++r)
            for (int t = 0; t < 3; ++t)
                res.WX[q].tail(r, t) = Tp[((q * npop + last) * R + r) * 3 + t];
    }
    return res;
}

/**
 * Schweitzer fixed point counterpart of pfqn_mvasjn.
 *
 * The closure is applied to the SIZE-RESOLVED queue length lam_k W_k(x) f_k(x)
 * rather than to its integral: removing one customer of class r scales the
 * class-r density by (N_r-1)/N_r and leaves the other classes unchanged.
 * Integrating over x recovers the usual Schweitzer rule, so the closure is the
 * exact analogue of the one used at the ordinary stations. Cost per iteration
 * is O(M R ns) against the prod(N+1) M R ns of the lattice, and the population
 * may be arbitrarily large. What is given up is the population dependence of
 * the SHAPE of W(x): its level may scale but its shape is fixed, whereas the
 * true profile stiffens with the load. The error therefore concentrates at high
 * utilization, where the SJN approximation is already weakest.
 */
inline SjnResult pfqn_amvasjn(const Matrix<double>& L, const std::vector<double>& N,
                              const std::vector<double>& Z, const Matrix<double>& scv,
                              const std::vector<std::size_t>& sjnset, const Matrix<double>& V,
                              const SjnOptions& options) {
    const detail::SjnArgs a = detail::sjn_args(L, N, Z, scv, sjnset, V, options);
    const std::size_t M = a.M, R = a.R, nsjn = a.sjnset.size();
    const bool useprio = !a.opt.prio.empty();
    const std::size_t ns = a.opt.ns, ngrid = ns + 1;

    std::vector<detail::SjnGrid> G(nsjn);
    for (std::size_t q = 0; q < nsjn; ++q) {
        std::vector<double> Sq(R), scvq(R);
        for (std::size_t r = 0; r < R; ++r) {
            Sq[r] = a.S(a.sjnset[q], r);
            scvq[r] = a.scv(a.sjnset[q], r);
        }
        G[q] = detail::sjn_setup(Sq, scvq, ns, a.opt.Lfactor);
    }

    // start from the product-form Schweitzer solution: a light-load guess would
    // put the deflated utilization above one, where the SJN denominator has no
    // solution at all
    const AmvaResult<double> bs =
        pfqn_bs(L, a.N, a.Z, std::vector<AmvaSched>(), a.opt.tol, a.opt.iter_max);
    SjnResult res;
    res.XN = bs.XN;
    res.QN = bs.QN;
    res.UN = bs.UN;
    res.CN = bs.RN;

    std::vector<Matrix<double>> W(nsjn, Matrix<double>(ngrid, R, 0.0));
    std::vector<Matrix<double>> P(nsjn, Matrix<double>(ngrid, R, 0.0));
    Matrix<double> Iinf(nsjn, R, 0.0);
    std::vector<Matrix<double>> Tail(nsjn, Matrix<double>(R, 3, 0.0));

    std::size_t it = 0;
    bool converged = false;
    double delta = 0.0;
    while (!converged && it < a.opt.iter_max) {
        ++it;
        Matrix<double> Cit = res.CN;
        std::vector<Matrix<double>> Wit = W, Pit = P, Tit = Tail;
        Matrix<double> Iit = Iinf;
        for (std::size_t r = 0; r < R; ++r) {
            if (a.N[r] == 0.0) continue;
            std::vector<double> beta(R, 1.0);
            beta[r] = (a.N[r] - 1.0) / a.N[r];
            for (std::size_t m = 0; m < M; ++m) {
                std::size_t q = nsjn;
                for (std::size_t t = 0; t < nsjn; ++t)
                    if (a.sjnset[t] == m) q = t;
                if (q == nsjn) {
                    double qsum = 0.0;
                    for (std::size_t k = 0; k < R; ++k) qsum += beta[k] * res.QN(m, k);
                    Cit(m, r) = L(m, r) * (1.0 + qsum);
                    continue;
                }
                detail::SjnState st;
                st.lam.assign(R, 0.0);
                st.U.assign(R, 0.0);
                st.Q.assign(R, 0.0);
                st.phiinf.assign(R, 0.0);
                for (std::size_t k = 0; k < R; ++k) {
                    st.lam[k] = res.XN[k] * a.V(m, k);
                    st.U[k] = res.UN(m, k);
                    st.Q[k] = res.QN(m, k);
                    st.phiinf[k] = Iinf(q, k);
                }
                st.W = &W[q];
                st.phi = &P[q];
                std::vector<double> Sm(R), scvm(R), Vm(R);
                for (std::size_t k = 0; k < R; ++k) {
                    Sm[k] = a.S(m, k);
                    scvm[k] = a.scv(m, k);
                    Vm[k] = a.V(m, k);
                }
                const detail::SjnStationResult sr =
                    detail::sjn_station(m, r, G[q], Sm, scvm, Vm, st, beta, useprio, a.opt.prio);
                Cit(m, r) = sr.C;
                for (std::size_t i = 0; i < ngrid; ++i) {
                    Wit[q](i, r) = sr.W[i];
                    Pit[q](i, r) = sr.phi[i];
                }
                Iit(q, r) = sr.phiinf;
                for (int t = 0; t < 3; ++t) Tit[q](r, t) = sr.tail[t];
            }
        }
        const detail::SjnCapResult cap = detail::sjn_cap(Cit, L, a.N, a.Z, a.sjnset, a.opt.umax);
        if (cap.bound) {
            res.capped = true;
            for (std::size_t q = 0; q < nsjn; ++q) {
                const double kq = cap.kappa[q];
                for (std::size_t i = 0; i < ngrid; ++i)
                    for (std::size_t r = 0; r < R; ++r) {
                        Wit[q](i, r) *= kq;
                        Pit[q](i, r) *= kq;
                    }
                for (std::size_t r = 0; r < R; ++r) {
                    Iit(q, r) *= kq;
                    Tit[q](r, 0) *= kq;
                    Tit[q](r, 1) *= kq;
                }
            }
        }
        Matrix<double> Qit(M, R, 0.0), Uit(M, R, 0.0);
        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t r = 0; r < R; ++r) {
                Qit(m, r) = cap.X[r] * Cit(m, r);
                Uit(m, r) = cap.X[r] * L(m, r);
            }
        delta = 0.0;
        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t r = 0; r < R; ++r)
                delta = std::max(delta, std::fabs(Qit(m, r) - res.QN(m, r)));
        for (std::size_t q = 0; q < nsjn; ++q)
            for (std::size_t i = 0; i < ngrid; ++i)
                for (std::size_t r = 0; r < R; ++r)
                    delta = std::max(delta, std::fabs(Wit[q](i, r) - W[q](i, r)));
        res.XN = cap.X;
        res.QN = Qit;
        res.UN = Uit;
        res.CN = Cit;
        W = Wit;
        P = Pit;
        Iinf = Iit;
        Tail = Tit;
        converged = delta < a.opt.tol;
    }
    res.iter = it;
    res.converged = converged;

    res.WX.resize(nsjn);
    for (std::size_t q = 0; q < nsjn; ++q) {
        res.WX[q].station = a.sjnset[q] + 1;
        res.WX[q].x = G[q].x;
        res.WX[q].W = W[q];
        res.WX[q].tail = Tail[q];
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SJN_H
