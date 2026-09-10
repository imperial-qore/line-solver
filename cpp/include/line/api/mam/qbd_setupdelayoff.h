/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_QBD_SETUPDELAYOFF_H
#define LINE_API_MAM_QBD_SETUPDELAYOFF_H

/**
 * Mean queue length of an M/M/1 queue with a setup delay and a delayed-off
 * period, solved as a QBD.
 *
 * Templated port of matlab/src/api/mam/qbd_setupdelayoff.m. The server is
 * switched off when the system empties, but only after a delay-off period of
 * rate betarate and SCV betascv has elapsed; an arrival during that period
 * finds the server still up. Once off, an arrival starts a setup of rate
 * alpharate and SCV alphascv before service can begin.
 *
 * The QBD phase index is overloaded by level, which is why the phases cannot
 * be redistributed on a level change:
 *   - at level 0 phase 1 is "server off" and phases na+1..na+nb are the
 *     delay-off phases;
 *   - above level 0 phases 1..na are the setup phases and phase na+1 is the
 *     busy server.
 * An arrival to an off server must therefore enter the setup at phase 1,
 * which is why both phases are built in CANONICAL COXIAN form: its entry
 * vector is [1 0 ... 0] for every SCV. The reference notes that
 * APH.fitMeanAndSCV violates this for SCV > 1, returning a hyperexponential
 * entered at phase 2 with probability 3/4, which the chain then silently
 * entered at phase 1.
 *
 * The phases are given as RATES and an exponential phase is built from its
 * rate directly rather than round-tripped through its mean, which is what the
 * reference does after the round trip turned a finite 1e8 rate into an
 * infinite one.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: the Coxian fit takes
 * a square root and R comes from cyclic reduction.
 *
 * TRUNCATION. The level series is cut where MATLAB's QBD_pi cuts it, at
 * accumulated mass 1 - 1e-10 or 501 level vectors, whichever comes first
 * (QBD_pi's MaxNumComp default is 500). qbd_pi's own default of 20000 levels
 * would keep more of the tail and report a slightly larger queue length, so
 * the cap is passed explicitly.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/qbd_r.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Sub-generator of the canonical Coxian form with the given RATE and SCV,
 * entered at phase 1 (the coxian_phase local of qbd_setupdelayoff.m, which
 * routes through Coxian.fitMeanAndSCV for SCV != 1).
 *
 * The four branches of Coxian.fitMeanAndSCV, with MEAN = 1/rate and
 * CoarseTol = 1e-3:
 *   |SCV - 1| <= tol      one exponential phase of rate 1/MEAN;
 *   0.5 + tol < SCV < 1 - tol
 *                         two phases, mu_i = 2/MEAN/(1 -+ sqrt(2 SCV - 1)),
 *                         phi = [0, 1], i.e. a pure series (hypoexponential);
 *   SCV <= 0.5 + tol      an Erlang of n = ceil(1/SCV) phases of rate n/MEAN;
 *   SCV > 1 + tol         two phases with mu_1 = 2/MEAN, mu_2 = mu_1/(2 SCV)
 *                         and phi_1 = 1 - mu_2/mu_1, the Coxian form of a
 *                         hyperexponential.
 * The sub-generator is diag(-mu) + diag(mu_i (1 - phi_i), 1), as in
 * Coxian.m's process assembly.
 */
template <class T>
Matrix<T> coxian_phase_subgen(const T& rate, const T& scv) {
    static_assert(num_traits<T>::has_transcendental,
                  "coxian_phase_subgen requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T zero = num_traits<T>::from_int(0);
    if (rate <= zero) throw InputError("coxian_phase_subgen: the rate must be positive");
    if (scv <= zero) throw InputError("coxian_phase_subgen: the SCV must be positive");

    // An exponential phase is built from the rate directly.
    if (scv == one) {
        Matrix<T> D0(1, 1);
        D0(0, 0) = -rate;
        return D0;
    }

    const T tol = num_traits<T>::from_double(1e-3);
    const T mean = one / rate;
    std::vector<T> mu, phi;
    if (scv >= one - tol && scv <= one + tol) {
        mu.push_back(one / mean);
        phi.push_back(one);
    } else if (scv > T(one / two) + tol && scv < one - tol) {
        using std::sqrt;
        const T s = T(sqrt(T(one + two * (scv - one))));
        mu.push_back(two / mean / (one + s));
        mu.push_back(two / mean / (one - s));
        phi.push_back(zero);
        phi.push_back(one);
    } else if (scv <= T(one / two) + tol) {
        const double inv = 1.0 / num_traits<T>::to_double(scv);
        const long n = static_cast<long>(std::ceil(inv));
        const T lambda = num_traits<T>::from_int(n) / mean;
        for (long k = 0; k < n; ++k) {
            mu.push_back(lambda);
            phi.push_back(zero);
        }
        phi[static_cast<std::size_t>(n) - 1] = one;
    } else {
        const T mu1 = two / mean;
        const T mu2 = mu1 / (two * scv);
        mu.push_back(mu1);
        mu.push_back(mu2);
        phi.push_back(T(one - mu2 / mu1));
        phi.push_back(one);
    }
    phi.back() = one;

    const std::size_t n = mu.size();
    Matrix<T> D0(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        D0(i, i) = -mu[i];
        if (i + 1 < n) D0(i, i + 1) = mu[i] * (one - phi[i]);
    }
    return D0;
}

/**
 * Mean queue length of the M/M/1 queue with setup delay and delay-off
 * (qbd_setupdelayoff.m).
 *
 * @param lambda    arrival rate
 * @param mu        service rate
 * @param alpharate rate of the setup phase
 * @param alphascv  SCV of the setup phase
 * @param betarate  rate of the delay-off phase
 * @param betascv   SCV of the delay-off phase
 */
template <class T>
T qbd_setupdelayoff(const T& lambda, const T& mu, const T& alpharate, const T& alphascv,
                    const T& betarate, const T& betascv) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_setupdelayoff requires transcendental arithmetic");
    using namespace qbd_detail;
    const T zero = num_traits<T>::from_int(0);

    const Matrix<T> Ta = coxian_phase_subgen(alpharate, alphascv);
    const std::size_t na = Ta.rows();
    std::vector<T> ta(na, zero);  // completion rate out of each setup phase
    for (std::size_t i = 0; i < na; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < na; ++j) s += Ta(i, j);
        ta[i] = -s;
    }

    const Matrix<T> Tb = coxian_phase_subgen(betarate, betascv);
    const std::size_t nb = Tb.rows();
    std::vector<T> tb(nb, zero);
    for (std::size_t i = 0; i < nb; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < nb; ++j) s += Tb(i, j);
        tb[i] = -s;
    }

    const std::size_t n = na + nb;
    Matrix<T> F(n, n, zero), B(n, n, zero), L(n, n, zero), L0(n, n, zero);

    for (std::size_t i = 0; i < na; ++i) F(i, i) = lambda;
    for (std::size_t i = 0; i < nb; ++i) F(na + i, na) = lambda;
    F(na, na) = lambda;
    B(na, na) = mu;

    for (std::size_t i = 0; i < na; ++i) {
        // whole generator row rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        for (std::size_t j = 0; j < na; ++j) L(i, j) = Ta(i, j);
        L(i, i) -= lambda;
        L(i, na) = ta[i];  // setup completes from phase i -> busy server
    }
    L(na, na) = -mu - lambda;
    for (std::size_t i = 1; i < nb; ++i) L(na + i, na + i) = -lambda;

    for (std::size_t i = 0; i < na; ++i) L0(i, i) = -lambda;
    for (std::size_t i = 0; i < nb; ++i) {
        for (std::size_t j = 0; j < nb; ++j) L0(na + i, na + j) = Tb(i, j);
        L0(na + i, na + i) -= lambda;
        L0(na + i, 0) = tb[i];  // delay-off expires from phase i -> server off
    }

    const QbdFundMat<T> fm = qbd_fundmat(B, L, F);
    const Matrix<T> pn =
        qbd_pi(B, L0, fm.R, static_cast<std::size_t>(501), T(num_traits<T>::from_double(1e-10)));

    // prior reference bug (15% high queue length): see _kb/03-api-layer.md (cpp port notes: mam)
    T QN = zero;
    for (std::size_t k = 1; k < pn.rows(); ++k) {
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += pn(k, j);
        QN += num_traits<T>::from_int(static_cast<long>(k)) * s;
    }
    return QN;
}

/** Mean queue length and throughput of the CLOSED setup/delay-off queue. */
template <class T>
struct SetupDelayoffClosed {
    T QN;  ///< mean number of jobs at the station
    T XN;  ///< throughput of the station
};

/**
 * Mean queue length and throughput of a FINITE-POPULATION queue with setup
 * delay and delay-off, a port of matlab/src/api/mam/qbd_setupdelayoff_closed.m.
 *
 * The closed twin of `qbd_setupdelayoff`. The population N is finite and Z is
 * the complementary delay, the mean time a customer spends away from this
 * station, so the arrival rate is state dependent, lambda(n) = (N - n)/Z, and
 * the level index is bounded by N. That makes the chain a LEVEL-DEPENDENT QBD
 * over finitely many levels, i.e. a finite CTMC, and it is solved exactly rather
 * than by a matrix-geometric tail.
 *
 * THE SEMANTICS ARE THE SIMULATOR'S, not the mean-value shortcut's. When the
 * queue empties the server begins a delay-off period; an arrival DURING it finds
 * the server still warm and resumes without setup (Solver_ssj's cancelDelayoff),
 * and only an arrival after the delay-off has expired pays the setup. That is an
 * M/M/1 with setup time AND close-down time. The per-instance cold-start race
 * `p_cold*E[setup] + S` this replaces raced the delay-off against the
 * per-instance idle time and carried NO queueing term, so it described a
 * serverless instance pool rather than a single-server vacation queue and left
 * the reported response time byte-identical across a tenfold change in the setup
 * mean.
 *
 * The phase index is overloaded by level exactly as in the open twin: at level 0
 * phase 1 is the OFF server and the rest are the delay-off; above level 0 the
 * phases are the setup and the last one is the busy server. Only the REACHABLE
 * states are enumerated, because a finite chain cannot carry an unreachable row:
 * it would be absorbing and the stationary solve singular.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental, like the open twin:
 * the Coxian fit takes a square root.
 *
 * @param N population of the closed chain
 * @param Z complementary delay, the mean time a customer spends away
 * @param mu service rate of the station
 * @param alpharate rate of the setup phase
 * @param alphascv SCV of the setup phase
 * @param betarate rate of the delay-off phase
 * @param betascv SCV of the delay-off phase
 */
template <class T>
SetupDelayoffClosed<T> qbd_setupdelayoff_closed(const T& N, const T& Z, const T& mu,
                                                const T& alpharate, const T& alphascv,
                                                const T& betarate, const T& betascv) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_setupdelayoff_closed requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    SetupDelayoffClosed<T> out;
    out.QN = zero;
    out.XN = zero;
    const long pop = static_cast<long>(std::llround(num_traits<T>::to_double(N)));
    if (pop <= 0 || num_traits<T>::to_double(mu) <= 0) return out;
    T Zc = Z;
    if (num_traits<T>::to_double(Zc) < lang::GlobalConstants::FineTol)
        Zc = num_traits<T>::from_double(lang::GlobalConstants::FineTol);

    const Matrix<T> Ta = coxian_phase_subgen(alpharate, alphascv);
    const std::size_t na = Ta.rows();
    std::vector<T> ta(na, zero);
    for (std::size_t i = 0; i < na; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < na; ++j) s += Ta(i, j);
        ta[i] = T(-s);
    }
    const Matrix<T> Tb = coxian_phase_subgen(betarate, betascv);
    const std::size_t nb = Tb.rows();
    std::vector<T> tb(nb, zero);
    for (std::size_t i = 0; i < nb; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < nb; ++j) s += Tb(i, j);
        tb[i] = T(-s);
    }

    const std::size_t off = 0;      // level 0, server off
    const std::size_t base = 1 + nb;  // level 0 delay-off occupies 1..nb
    const std::size_t P = na + 1;
    const std::size_t m = base + static_cast<std::size_t>(pop) * P;
    Matrix<T> Q(m, m, zero);

    // lambda(n) = (N - n)/Z, zero at the top level
    std::vector<T> lam(static_cast<std::size_t>(pop) + 1, zero);
    for (long n = 0; n < pop; ++n)
        lam[static_cast<std::size_t>(n)] = T(num_traits<T>::from_int(pop - n) / Zc);

    if (num_traits<T>::to_double(lam[0]) > 0) Q(off, base) += lam[0];
    for (std::size_t j = 0; j < nb; ++j) {
        const std::size_t rj = 1 + j;
        for (std::size_t j2 = 0; j2 < nb; ++j2)
            if (j2 != j) Q(rj, 1 + j2) += Tb(j, j2);
        Q(rj, off) += tb[j];
        // an arrival during the delay-off cancels it and resumes WITHOUT setup
        if (num_traits<T>::to_double(lam[0]) > 0) Q(rj, base + na) += lam[0];
    }
    for (long n = 1; n <= pop; ++n) {
        const std::size_t lvl = base + static_cast<std::size_t>(n - 1) * P;
        const std::size_t up = base + static_cast<std::size_t>(n) * P;
        const T& ln = lam[static_cast<std::size_t>(n)];
        const bool rising = n < pop && num_traits<T>::to_double(ln) > 0;
        for (std::size_t i = 0; i < na; ++i) {
            for (std::size_t i2 = 0; i2 < na; ++i2)
                if (i2 != i) Q(lvl + i, lvl + i2) += Ta(i, i2);
            Q(lvl + i, lvl + na) += ta[i];
            // an arrival during the setup joins the queue and the setup carries on
            // in the SAME phase: the level rises, the phase does not move
            if (rising) Q(lvl + i, up + i) += ln;
        }
        if (rising) Q(lvl + na, up + na) += ln;
        // a completion that empties the queue starts the delay-off at its phase 1
        const std::size_t down = n - 1 >= 1 ? base + static_cast<std::size_t>(n - 2) * P + na : 1;
        Q(lvl + na, down) += mu;
    }
    for (std::size_t i = 0; i < m; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < m; ++j)
            if (j != i) s += Q(i, j);
        Q(i, i) = T(-s);
    }

    const std::vector<T> pi = mc::ctmc_solve(Q);
    double total = 0.0;
    for (std::size_t i = 0; i < pi.size(); ++i)
        total += std::max(0.0, num_traits<T>::to_double(pi[i]));
    if (total <= 0) return out;

    T QN = zero, pbusy = zero;
    for (long n = 1; n <= pop; ++n) {
        const std::size_t lvl = base + static_cast<std::size_t>(n - 1) * P;
        T level = pi[lvl + na];
        for (std::size_t i = 0; i < na; ++i) level += pi[lvl + i];
        QN += T(num_traits<T>::from_int(n) * level);
        pbusy += pi[lvl + na];
    }
    const T norm = num_traits<T>::from_double(total);
    out.QN = T(QN / norm);
    out.XN = T(mu * pbusy / norm);
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_QBD_SETUPDELAYOFF_H
