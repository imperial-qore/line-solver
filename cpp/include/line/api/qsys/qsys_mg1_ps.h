/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1_PS_H
#define LINE_API_QSYS_QSYS_MG1_PS_H

/**
 * Sojourn-time distribution of the M/G/1 processor-sharing queue.
 *
 * Port of matlab/src/api/qsys/qsys_mg1_ps.m (twins `Qsys_mg1_ps.java`,
 * `python/line_solver/api/qsys/mg1ps.py`). MATLAB is the reference.
 *
 * Jobs arrive Poisson at rate lambda at one egalitarian processor-sharing
 * server whose service requirement has LST bhat(tau) and mean m1. Writing V(x)
 * for the sojourn of a tagged job of requirement x and rho = lambda m1 < 1, Ott
 * (1984) and Yashkov (1983) give
 *
 *     E[exp(-s V(x))] = (1 - rho) / D(s,x),
 *
 * D(s,x) being the inverse Laplace transform, evaluated at x, of
 *
 *     f(tau;s) = [ (1-rho) tau^2 - (1-rho) lambda (1-bhat(tau)) tau
 *                  + s rho tau - s lambda (1-bhat(tau)) ]
 *                / [ tau^2 (tau - s - lambda (1-bhat(tau))) ].
 *
 * The transform is exact but IMPLICIT: f has to be inverted in tau. For
 * phase-type service f is a proper rational function of tau, the double pole at
 * the origin cancels, and D(s,x) comes out in closed form as a finite sum of
 * residues -- so the M/PH/1-PS queue, hence every service law that can be
 * fitted phase-type, is exactly solvable. For a transform supplied as a
 * callback, f is inverted numerically on a Bromwich contour placed to the right
 * of the dominant singularity tau*(s), the unique root of
 * tau = s + lambda (1 - bhat(tau)) in the right half plane, which the fixed
 * point of that equation reaches at geometric rate rho.
 *
 * THE CONDITIONAL SOJOURN IS ATOMIC ON THE LATTICE t = (k+1) x, and that is not
 * a numerical artifact. Processor sharing gives every job in the system the
 * same service rate, so if the k jobs present on arrival all outlive the tagged
 * job and nothing arrives, the sojourn is exactly (k+1) x. The k = 0 atom,
 * (1-rho) exp(-lambda x), is the probability of finding the system empty and
 * sharing it with nobody, and it is the only one that stays exact for general
 * service. It is removed before inverting in s; the remaining atoms make the
 * conditional CDF JUMP and leave no density, so the density is reported as NaN
 * on the lattice rather than as the finite garbage a smooth inversion returns
 * there.
 *
 * TWO SUBSTITUTIONS FOR MATLAB, both self-contained rather than approximate:
 *
 * 1. `roots(P)` is computed by Durand-Kerner rather than by the eigenvalues of
 *    the companion matrix, so this header does not require LAPACK. The
 *    polynomial has degree n+1 in the phase count, which is small, and the
 *    iteration converges to machine precision on it; the roots are a SET, and
 *    the residue sum that consumes them does not depend on their order.
 * 2. The Golub-Welsch call behind the panelled quadrature is the tree's own
 *    `gauss_legendre`, which computes the same Legendre nodes by Newton. The
 *    reference notes that the rule is fixed so that every codebase shares it;
 *    the nodes agree to machine precision, so they do.
 *
 * ARITHMETIC: double. The inversion, the root finding and the contour are all
 * inherently floating point, and the reference's own tolerances are absolute in
 * double.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

using Cplx = std::complex<double>;

/** Options of `qsys_mg1_ps`, mirroring the reference's name-value pairs. */
struct Mg1PsOptions {
    std::vector<double> x;   ///< service requirements to condition on
    std::vector<double> s;   ///< transform arguments to tabulate the LST at
    std::vector<double> t;   ///< times to evaluate the sojourn distribution at
    std::size_t nterms = 41; ///< function evaluations per numerical inversion, ODD
    /** Service density, needed to remove the conditioning on a callback path. */
    std::function<double(double)> pdf;
};

/** Everything `qsys_mg1_ps` returns. */
struct Mg1PsResult {
    double rho = 0.0;
    double m1 = 0.0;
    double m2 = std::numeric_limits<double>::quiet_NaN();
    std::function<Cplx(Cplx, double)> lstCond;    ///< (s,x) -> E[exp(-s V(x))]
    std::function<Cplx(Cplx, double)> lstExcess;  ///< (s,x) -> E[exp(-s (V(x)-x))]
    std::function<Cplx(Cplx)> lstUncond;          ///< s -> E[exp(-s V)]
    std::function<Cplx(Cplx)> dominantRoot;       ///< s -> tau*(s)
    std::vector<double> x, s, t;
    Matrix<double> lstCondVal;                    ///< (|x| x |s|)
    std::vector<double> lstUncondVal;
    std::vector<double> atomCond;                 ///< (1-rho) exp(-lambda x), at t = x
    double atomUncond = 0.0;                      ///< (1-rho) bhat(lambda)
    std::vector<double> meanCond, m2Cond, varCond;
    double meanUncond = 0.0;
    double m2Uncond = std::numeric_limits<double>::quiet_NaN();
    double varUncond = std::numeric_limits<double>::quiet_NaN();
    Matrix<double> pdfCond, cdfCond;              ///< (|x| x |t|)
    std::vector<double> pdfUncond, cdfUncond;
};

namespace mg1psdetail {

/** Polynomial product, highest degree first, as MATLAB's `conv`. */
inline std::vector<double> conv(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.empty() || b.empty()) return std::vector<double>();
    std::vector<double> c(a.size() + b.size() - 1, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
    return c;
}

/** Horner evaluation of a polynomial given highest degree first. */
template <class C, class U>
U polyval(const std::vector<C>& p, const U& z) {
    U v = U(0.0);
    for (std::size_t i = 0; i < p.size(); ++i) v = v * z + U(p[i]);
    return v;
}

/** Derivative of a polynomial given highest degree first. */
inline std::vector<double> polyder(const std::vector<double>& p) {
    if (p.size() <= 1) return std::vector<double>(1, 0.0);
    const std::size_t n = p.size() - 1;
    std::vector<double> d(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) d[i] = p[i] * static_cast<double>(n - i);
    return d;
}

/**
 * Roots of a polynomial by Durand-Kerner.
 *
 * MATLAB uses the companion-matrix eigenvalues; this is LAPACK-free and exact
 * to machine precision at the degrees involved (one more than the phase count).
 * The roots are a SET and every consumer here is symmetric in them.
 */
inline std::vector<Cplx> poly_roots(const std::vector<double>& p_in) {
    std::vector<double> p = p_in;
    // Strip leading zeros: they are not roots, they are a lower degree.
    std::size_t lead = 0;
    while (lead + 1 < p.size() && p[lead] == 0.0) ++lead;
    p.erase(p.begin(), p.begin() + static_cast<long>(lead));
    if (p.size() <= 1) return std::vector<Cplx>();
    const std::size_t n = p.size() - 1;
    std::vector<double> mon(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) mon[i] = p[i] / p[0];

    // The classic spiral start, off the real axis so a real polynomial's
    // iterates do not collapse onto it.
    std::vector<Cplx> r(n);
    const Cplx seed(0.4, 0.9);
    Cplx pw(1.0, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        r[i] = pw;
        pw *= seed;
    }
    for (int it = 0; it < 1000; ++it) {
        double move = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            Cplx den(1.0, 0.0);
            for (std::size_t j = 0; j < n; ++j)
                if (j != i) den *= (r[i] - r[j]);
            if (std::abs(den) < 1e-300) continue;
            const Cplx step = polyval(mon, r[i]) / den;
            r[i] -= step;
            move = std::max(move, std::abs(step));
        }
        if (move < 1e-14) break;
    }
    return r;
}

/**
 * Matrix exponential of a COMPLEX matrix, by scaling and squaring with a
 * Taylor series.
 *
 * `util/expm.h` is templated on `num_traits`, which has no complex
 * instantiation, so the repeated-pole branch below needs its own. Scaling to a
 * norm under 1/2 makes the truncated series converge to machine precision in
 * around twenty terms, and the squaring recovers the full argument.
 */
inline Matrix<Cplx> expm_cplx(const Matrix<Cplx>& A) {
    const std::size_t n = A.rows();
    double nrm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double row = 0.0;
        for (std::size_t j = 0; j < n; ++j) row += std::abs(A(i, j));
        nrm = std::max(nrm, row);
    }
    int sq = 0;
    while (nrm > 0.5) {
        nrm /= 2.0;
        ++sq;
    }
    const double sc = std::pow(2.0, -static_cast<double>(sq));
    Matrix<Cplx> B(n, n, Cplx(0.0, 0.0)), E(n, n, Cplx(0.0, 0.0)), Tk(n, n, Cplx(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        E(i, i) = Cplx(1.0, 0.0);
        Tk(i, i) = Cplx(1.0, 0.0);
        for (std::size_t j = 0; j < n; ++j) B(i, j) = A(i, j) * sc;
    }
    for (int k = 1; k <= 30; ++k) {
        Matrix<Cplx> P(n, n, Cplx(0.0, 0.0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                Cplx acc(0.0, 0.0);
                for (std::size_t l = 0; l < n; ++l) acc += Tk(i, l) * B(l, j);
                P(i, j) = acc / static_cast<double>(k);
            }
        Tk = P;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) E(i, j) += Tk(i, j);
    }
    for (int t = 0; t < sq; ++t) {
        Matrix<Cplx> P(n, n, Cplx(0.0, 0.0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                Cplx acc(0.0, 0.0);
                for (std::size_t l = 0; l < n; ++l) acc += E(i, l) * E(l, j);
                P(i, j) = acc;
            }
        E = P;
    }
    return E;
}

/**
 * Abate-Whitt Euler inversion, symmetrized so that a complex-valued time
 * function is handled as well as a real-valued one.
 */
template <class Fn>
Cplx ilt(Fn fun, double t, std::size_t nterms) {
    const long ne = static_cast<long>((nterms - 1) / 2);
    std::vector<double> eta(static_cast<std::size_t>(2 * ne + 1), 0.0);
    eta[0] = 0.5;
    for (long i = 1; i <= ne; ++i) eta[static_cast<std::size_t>(i)] = 1.0;
    eta[static_cast<std::size_t>(2 * ne)] = std::pow(2.0, -static_cast<double>(ne));
    for (long k = 1; k <= ne - 1; ++k)
        eta[static_cast<std::size_t>(2 * ne - k)] =
            eta[static_cast<std::size_t>(2 * ne - k + 1)] +
            std::exp(std::lgamma(static_cast<double>(ne) + 1.0) -
                     static_cast<double>(ne) * std::log(2.0) -
                     std::lgamma(static_cast<double>(k) + 1.0) -
                     std::lgamma(static_cast<double>(ne - k) + 1.0));
    Cplx g(0.0, 0.0);
    const double pref = std::pow(10.0, static_cast<double>(ne) / 3.0);
    for (long k = 0; k <= 2 * ne; ++k) {
        const Cplx beta(static_cast<double>(ne) * std::log(10.0) / 3.0,
                        M_PI * static_cast<double>(k));
        const double sgn = (k % 2 == 0) ? 1.0 : -1.0;
        const double e = pref * sgn * eta[static_cast<std::size_t>(k)];
        const Cplx bj = beta / t;
        g += 0.5 * e * (fun(bj) + fun(std::conj(bj)));
    }
    return g / t;
}

/**
 * The panelled Gauss-Legendre rule the reference uses to remove the
 * conditioning: eight panels growing geometrically toward ymax, so that both
 * ends of an exponentially decaying density are resolved, times a 32-point rule.
 * The rule is FIXED, so it is the same in every codebase.
 */
inline void quad_nodes(double ymax, std::vector<double>* y, std::vector<double>* w,
                       std::size_t npanel = 8, std::size_t ng = 32) {
    y->clear();
    w->clear();
    std::vector<double> edges;
    edges.push_back(0.0);
    for (long k = -static_cast<long>(npanel); k <= 0; ++k)
        edges.push_back(ymax * std::pow(2.0, static_cast<double>(k)));
    for (std::size_t k = 0; k + 1 < edges.size(); ++k) {
        std::vector<double> xg, wg;
        pfqn::detail::gauss_legendre<double>(ng, edges[k], edges[k + 1], xg, wg);
        for (std::size_t i = 0; i < xg.size(); ++i) {
            y->push_back(xg[i]);
            w->push_back(wg[i]);
        }
    }
}

/**
 * Second moment from the transform's curvature at the origin.
 *
 * The stencil is ONE-SIDED so the transform is never sampled at a negative
 * argument, where it need not converge; the step is scaled by the conditional
 * mean but capped by the unconditional one so it stays finite as x -> 0; and
 * Richardson extrapolation over h and h/2 removes the leading truncation.
 */
template <class Fn>
double second_moment(Fn lst, double meanref, double meanscale) {
    if (!(meanref > 0.0)) return 0.0;
    const double h = std::min(1e-2 / meanref, 1.0 / meanscale);
    auto d2 = [&lst](double hh) {
        double f[6];
        for (int j = 0; j < 6; ++j) f[j] = std::real(lst(Cplx(static_cast<double>(j) * hh, 0.0)));
        return (45.0 * f[0] - 154.0 * f[1] + 214.0 * f[2] - 156.0 * f[3] + 61.0 * f[4] -
                10.0 * f[5]) /
               (12.0 * hh * hh);
    };
    return (16.0 * d2(h / 2.0) - d2(h)) / 15.0;
}

}  // namespace mg1psdetail

/**
 * @param lambda arrival rate, finite and positive
 * @param alpha  phase-type initial probability vector
 * @param Tmat   phase-type subgenerator
 * @param opt    grids and tuning
 */
inline Mg1PsResult qsys_mg1_ps(double lambda, const std::vector<double>& alpha,
                               const Matrix<double>& Tmat, const Mg1PsOptions& opt) {
    using namespace mg1psdetail;
    if (!(lambda > 0.0) || !std::isfinite(lambda))
        throw InputError("qsys_mg1_ps: lambda must be a finite positive scalar");
    if (opt.nterms % 2 == 0 || opt.nterms < 11)
        throw InputError("qsys_mg1_ps: nterms must be an odd integer of at least 11");

    const std::size_t n = alpha.size();
    if (n == 0 || Tmat.rows() != n || Tmat.cols() != n)
        throw InputError("qsys_mg1_ps: T must be n x n to match alpha");
    double asum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (alpha[i] < -1e-12) throw InputError("qsys_mg1_ps: alpha must be a probability vector");
        asum += alpha[i];
    }
    if (std::fabs(asum - 1.0) > 1e-8)
        throw InputError("qsys_mg1_ps: alpha must be a probability vector");
    std::vector<double> exitrate(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double row = 0.0;
        for (std::size_t j = 0; j < n; ++j) row += Tmat(i, j);
        exitrate[i] = -row;
        if (exitrate[i] < -1e-10 || Tmat(i, i) >= 0.0)
            throw InputError("qsys_mg1_ps: T must be a proper phase-type subgenerator");
    }

    // Moments, by solving with (-T) rather than forming its inverse.
    Matrix<double> negT(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negT(i, j) = -Tmat(i, j);
    const std::vector<double> ones(n, 1.0);
    const std::vector<double> u1 = solve(negT, ones);
    const std::vector<double> u2 = solve(negT, u1);
    double m1 = 0.0, m2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        m1 += alpha[i] * u1[i];
        m2 += 2.0 * alpha[i] * u2[i];
    }
    const double rho = lambda * m1;
    if (rho >= 1.0)
        throw InputError("qsys_mg1_ps: the system is unstable, the utilization is at least one");

    // Faddeev-LeVerrier gives det(tau I - T) and the adjugate in one sweep, so
    // bhat(tau) = nb(tau)/db(tau) as polynomials of degree n-1 and n.
    std::vector<double> db(n + 1, 0.0), nb(n, 0.0);
    db[0] = 1.0;
    Matrix<double> Mk(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) Mk(i, i) = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        double v = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double row = 0.0;
            for (std::size_t j = 0; j < n; ++j) row += Mk(i, j) * exitrate[j];
            v += alpha[i] * row;
        }
        nb[k - 1] = v;
        Matrix<double> TM(n, n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                double acc = 0.0;
                for (std::size_t l = 0; l < n; ++l) acc += Tmat(i, l) * Mk(l, j);
                TM(i, j) = acc;
            }
        double tr = 0.0;
        for (std::size_t i = 0; i < n; ++i) tr += TM(i, i);
        db[k] = -tr / static_cast<double>(k);
        Mk = TM;
        for (std::size_t i = 0; i < n; ++i) Mk(i, i) += db[k];
    }

    auto bhat = [db, nb](const Cplx& tau) { return polyval(nb, tau) / polyval(db, tau); };
    auto bpdf = [&Tmat, &alpha, &exitrate, n](double y) {
        Matrix<double> E = Tmat;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) E(i, j) = Tmat(i, j) * y;
        const Matrix<double> Ey = expm(E);
        double v = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double row = 0.0;
            for (std::size_t j = 0; j < n; ++j) row += Ey(i, j) * exitrate[j];
            v += alpha[i] * row;
        }
        return v;
    };

    // D(s,x) = exp(scale x) val, exactly, from the residues of f(tau;s), whose
    // double pole at the origin cancels. The scale is carried separately so
    // that neither factor overflows at large s.
    auto denom = [db, nb, n, lambda, rho](const Cplx& s, double x, double* scale_out) {
        std::vector<double> dm(db.size(), 0.0);
        dm[0] = db[0];
        for (std::size_t i = 1; i < db.size(); ++i) dm[i] = db[i] - nb[i - 1];

        // These carry s, so they are complex; the polynomial helpers are
        // written over a scalar type for exactly that reason.
        std::vector<Cplx> Pc, Ac;
        {
            std::vector<Cplx> lin;
            lin.push_back(Cplx(1.0, 0.0));
            lin.push_back(-(s + Cplx(lambda, 0.0)));
            std::vector<Cplx> dbz(db.size());
            for (std::size_t i = 0; i < db.size(); ++i) dbz[i] = Cplx(db[i], 0.0);
            Pc.assign(lin.size() + dbz.size() - 1, Cplx(0.0, 0.0));
            for (std::size_t i = 0; i < lin.size(); ++i)
                for (std::size_t j = 0; j < dbz.size(); ++j) Pc[i + j] += lin[i] * dbz[j];
            // + [0, 0, lambda*nb]
            for (std::size_t i = 0; i < nb.size(); ++i)
                Pc[Pc.size() - nb.size() + i] += Cplx(lambda * nb[i], 0.0);

            std::vector<Cplx> q1;
            q1.push_back(Cplx(1.0 - rho, 0.0));
            q1.push_back(s * Cplx(rho, 0.0));
            q1.push_back(Cplx(0.0, 0.0));
            std::vector<Cplx> t1(q1.size() + dbz.size() - 1, Cplx(0.0, 0.0));
            for (std::size_t i = 0; i < q1.size(); ++i)
                for (std::size_t j = 0; j < dbz.size(); ++j) t1[i + j] += q1[i] * dbz[j];

            std::vector<Cplx> q2;
            q2.push_back(Cplx(1.0 - rho, 0.0));
            q2.push_back(s);
            std::vector<Cplx> ldm(dm.size());
            for (std::size_t i = 0; i < dm.size(); ++i) ldm[i] = Cplx(lambda * dm[i], 0.0);
            std::vector<Cplx> t2(q2.size() + ldm.size() - 1, Cplx(0.0, 0.0));
            for (std::size_t i = 0; i < q2.size(); ++i)
                for (std::size_t j = 0; j < ldm.size(); ++j) t2[i + j] += q2[i] * ldm[j];

            Ac.assign(t1.size(), Cplx(0.0, 0.0));
            for (std::size_t i = 0; i < t1.size(); ++i) Ac[i] = t1[i];
            for (std::size_t i = 0; i < t2.size(); ++i) Ac[Ac.size() - t2.size() + i] -= t2[i];
        }
        // The double pole at the origin must cancel: the two lowest
        // coefficients of A are what would survive it.
        double nrmA = 0.0;
        for (std::size_t i = 0; i < Ac.size(); ++i) nrmA = std::max(nrmA, std::abs(Ac[i]));
        const double tailA = std::abs(Ac[Ac.size() - 2]) + std::abs(Ac[Ac.size() - 1]);
        if (tailA > 1e-6 * std::max(1.0, nrmA))
            throw InputError("qsys_mg1_ps: the double pole at the origin did not cancel");
        // `A(1:n+1)` in the reference: the first n+1 of the n+3 coefficients.
        // The two dropped ones are the cancelled double pole, checked just
        // above; taking n+2 instead silently keeps half of it and the LST comes
        // back above one.
        std::vector<Cplx> Ahat(Ac.begin(), Ac.begin() + static_cast<long>(n + 1));

        // Durand-Kerner over the complex coefficients of P.
        std::vector<Cplx> r;
        {
            std::vector<Cplx> mon(Pc.size());
            for (std::size_t i = 0; i < Pc.size(); ++i) mon[i] = Pc[i] / Pc[0];
            const std::size_t deg = mon.size() - 1;
            r.assign(deg, Cplx(0.0, 0.0));
            const Cplx seed(0.4, 0.9);
            Cplx pw(1.0, 0.0);
            for (std::size_t i = 0; i < deg; ++i) {
                r[i] = pw;
                pw *= seed;
            }
            for (int it = 0; it < 2000; ++it) {
                double move = 0.0;
                for (std::size_t i = 0; i < deg; ++i) {
                    Cplx den(1.0, 0.0);
                    for (std::size_t j = 0; j < deg; ++j)
                        if (j != i) den *= (r[i] - r[j]);
                    if (std::abs(den) < 1e-300) continue;
                    const Cplx step = polyval(mon, r[i]) / den;
                    r[i] -= step;
                    move = std::max(move, std::abs(step));
                }
                if (move < 1e-14) break;
            }
        }
        double scale = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < r.size(); ++i) scale = std::max(scale, r[i].real());
        *scale_out = scale;

        double maxr = 0.0, minsep = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < r.size(); ++i) maxr = std::max(maxr, std::abs(r[i]));
        for (std::size_t i = 0; i < r.size(); ++i)
            for (std::size_t j = 0; j < r.size(); ++j)
                if (i != j) minsep = std::min(minsep, std::abs(r[i] - r[j]));

        if (r.size() >= 2 && minsep > 1e-7 * std::max(1.0, maxr)) {
            std::vector<Cplx> dP(Pc.size() - 1, Cplx(0.0, 0.0));
            for (std::size_t i = 0; i + 1 < Pc.size(); ++i)
                dP[i] = Pc[i] * static_cast<double>(Pc.size() - 1 - i);
            Cplx val(0.0, 0.0);
            for (std::size_t i = 0; i < r.size(); ++i) {
                const Cplx coef = polyval(Ahat, r[i]) / polyval(dP, r[i]);
                val += std::exp((r[i] - Cplx(scale, 0.0)) * x) * coef;
            }
            return val;
        }
        // Repeated poles: the companion realization of Ahat/P, exactly as the
        // reference falls back to. A residue sum is undefined there; the
        // matrix exponential is not.
        const std::size_t d = Pc.size() - 1;
        Matrix<Cplx> Acomp(d, d, Cplx(0.0, 0.0));
        for (std::size_t j = 0; j < d; ++j) Acomp(0, j) = -Pc[j + 1] / Pc[0];
        for (std::size_t i = 1; i < d; ++i) Acomp(i, i - 1) = Cplx(1.0, 0.0);
        for (std::size_t i = 0; i < d; ++i) Acomp(i, i) -= Cplx(scale, 0.0);
        Matrix<Cplx> Ax(d, d, Cplx(0.0, 0.0));
        for (std::size_t i = 0; i < d; ++i)
            for (std::size_t j = 0; j < d; ++j) Ax(i, j) = Acomp(i, j) * x;
        const Matrix<Cplx> E = expm_cplx(Ax);
        Cplx val(0.0, 0.0);
        for (std::size_t j = 0; j < d && j < Ahat.size(); ++j) val += (Ahat[j] / Pc[0]) * E(j, 0);
        return val;
    };

    Mg1PsResult res;
    res.rho = rho;
    res.m1 = m1;
    res.m2 = m2;
    res.meanUncond = m1 / (1.0 - rho);
    res.atomUncond = (1.0 - rho) * bhat(Cplx(lambda, 0.0)).real();

    auto lst_cond = [denom, rho](const Cplx& s, double x, bool excess) {
        if (x == 0.0) return Cplx(1.0, 0.0);
        double scale = 0.0;
        const Cplx val = denom(s, x, &scale);
        const Cplx expo = (excess ? s : Cplx(0.0, 0.0)) - Cplx(scale, 0.0);
        return Cplx(1.0 - rho, 0.0) * std::exp(expo * x) / val;
    };
    res.lstCond = [lst_cond](Cplx s, double x) { return lst_cond(s, x, false); };
    res.lstExcess = [lst_cond](Cplx s, double x) { return lst_cond(s, x, true); };
    res.dominantRoot = [lambda, bhat, rho](Cplx s) {
        Cplx tau = s;
        const int maxit = std::max(200, static_cast<int>(std::ceil(
                                            3.0 * std::log(1e-15) / std::log(std::max(rho, 1e-3)))));
        for (int it = 0; it < maxit; ++it) {
            const Cplx nx = s + Cplx(lambda, 0.0) * (Cplx(1.0, 0.0) - bhat(tau));
            if (std::abs(nx - tau) <= 1e-14 * std::max(1.0, std::abs(nx))) return nx;
            tau = nx;
        }
        throw InputError("qsys_mg1_ps: the dominant root iteration did not converge");
    };

    // The quadrature that removes the conditioning. The nodes do not move with
    // s, so the service density is sampled once.
    double mineig = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n; ++i) mineig = std::min(mineig, -Tmat(i, i));
    const double ymax = std::max(40.0 * m1, 40.0 / std::max(mineig, 1e-12));
    std::vector<double> yq, wq, bq;
    quad_nodes(ymax, &yq, &wq);
    bq.resize(yq.size());
    for (std::size_t k = 0; k < yq.size(); ++k) bq[k] = opt.pdf ? opt.pdf(yq[k]) : bpdf(yq[k]);
    res.lstUncond = [lst_cond, yq, wq, bq](Cplx s) {
        Cplx v(0.0, 0.0);
        for (std::size_t k = 0; k < yq.size(); ++k)
            v += wq[k] * bq[k] * lst_cond(s, yq[k], false);
        return v;
    };

    res.x = opt.x;
    res.s = opt.s;
    res.t = opt.t;
    res.lstCondVal = Matrix<double>(res.x.size(), res.s.size(), 0.0);
    for (std::size_t i = 0; i < res.x.size(); ++i)
        for (std::size_t j = 0; j < res.s.size(); ++j)
            res.lstCondVal(i, j) = res.lstCond(Cplx(res.s[j], 0.0), res.x[i]).real();
    res.lstUncondVal.assign(res.s.size(), 0.0);
    for (std::size_t j = 0; j < res.s.size(); ++j)
        res.lstUncondVal[j] = res.lstUncond(Cplx(res.s[j], 0.0)).real();

    res.atomCond.assign(res.x.size(), 0.0);
    res.meanCond.assign(res.x.size(), 0.0);
    res.m2Cond.assign(res.x.size(), 0.0);
    res.varCond.assign(res.x.size(), 0.0);
    for (std::size_t i = 0; i < res.x.size(); ++i) {
        res.atomCond[i] = (1.0 - rho) * std::exp(-lambda * res.x[i]);
        res.meanCond[i] = res.x[i] / (1.0 - rho);
        const double xi = res.x[i];
        res.m2Cond[i] = second_moment([&res, xi](Cplx u) { return res.lstCond(u, xi); },
                                      res.meanCond[i], res.meanUncond);
        res.varCond[i] = res.m2Cond[i] - res.meanCond[i] * res.meanCond[i];
    }

    // The unconditional second moment is integrated from the conditional one,
    // which is far better conditioned than differentiating the quadrature that
    // removed the conditioning.
    double m2u = 0.0;
    for (std::size_t k = 0; k < yq.size(); ++k) {
        const double yk = yq[k];
        m2u += wq[k] * bq[k] *
               second_moment([&res, yk](Cplx u) { return res.lstCond(u, yk); },
                             yk / (1.0 - rho), res.meanUncond);
    }
    res.m2Uncond = m2u;
    res.varUncond = m2u - res.meanUncond * res.meanUncond;

    res.pdfCond = Matrix<double>(res.x.size(), res.t.size(), 0.0);
    res.cdfCond = Matrix<double>(res.x.size(), res.t.size(), 0.0);
    for (std::size_t i = 0; i < res.x.size(); ++i) {
        const double atom = res.atomCond[i], xi = res.x[i];
        // V(x) >= x with an atom at x, so what is inverted is the EXCESS
        // V(x)-x net of its atom.
        auto gpdf = [&res, xi, atom](Cplx u) { return res.lstExcess(u, xi) - Cplx(atom, 0.0); };
        auto gcdf = [gpdf](Cplx u) { return gpdf(u) / u; };
        for (std::size_t j = 0; j < res.t.size(); ++j) {
            const double tj = res.t[j];
            if (tj < xi) continue;
            if (tj == xi) {
                res.cdfCond(i, j) = atom;
                continue;
            }
            res.cdfCond(i, j) = ilt(gcdf, tj - xi, opt.nterms).real() + atom;
            const double ratio = tj / xi;
            if (std::fabs(ratio - std::floor(ratio + 0.5)) < 1e-9)
                res.pdfCond(i, j) = std::numeric_limits<double>::quiet_NaN();
            else
                res.pdfCond(i, j) = ilt(gpdf, tj - xi, opt.nterms).real();
        }
    }

    res.pdfUncond.assign(res.t.size(), 0.0);
    res.cdfUncond.assign(res.t.size(), 0.0);
    for (std::size_t j = 0; j < res.t.size(); ++j) {
        if (res.t[j] <= 0.0) continue;
        res.pdfUncond[j] = ilt(res.lstUncond, res.t[j], opt.nterms).real();
        res.cdfUncond[j] =
            ilt([&res](Cplx u) { return res.lstUncond(u) / u; }, res.t[j], opt.nterms).real();
    }
    return res;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1_PS_H
