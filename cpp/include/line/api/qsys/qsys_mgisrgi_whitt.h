/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MGISRGI_WHITT_H
#define LINE_API_QSYS_MGISRGI_WHITT_H

/**
 * Engineering solution of the call-center model M/GI/s/r+GI.
 *
 * Templated port of matlab/src/api/qsys/qsys_mgisrgi_whitt.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_mgisrgi_whitt.java.
 *
 * Poisson arrivals at rate lambda, iid general service of mean 1/mu, s servers,
 * r extra waiting spaces and iid general patience.
 *
 * TWO APPROXIMATIONS. The general patience law becomes STATE-DEPENDENT
 * Markovian abandonment, a customer jth from the end of the queue abandoning at
 * rate delta_j = h(j/lambda) for the patience hazard h (eq. 3.3), because such a
 * customer has been waiting for about j/lambda; the total with k waiting is
 * Delta_k = sum_{j<=k} delta_j (eq. 3.4). The general service law becomes an
 * exponential of the same mean (Sec. 5). What is left is M/M/s/r+M(n), a
 * birth-and-death process:
 *
 *   mu_k = k mu                     1 <= k <= s
 *        = s mu + Delta_{k-s}       s+1 <= k <= s+r          (7.1)
 *   x_s = 1, x_{s+k+1} = lambda x_{s+k}/mu_{s+k+1}, x_{k-1} = mu_k x_k/lambda
 *   p_k = x_k / sum_j x_j                                    (7.4)-(7.7)
 *
 * and the customer experience follows from the kernel
 *
 *   m_k(j)   = 1/(s mu + Delta_k - Delta_{j-1})              (7.11)
 *   phi_k(j) = delta_j m_k(j)                                (7.10)
 *   sigma_k  = prod_j (1 - phi_k(j))                         (7.13)
 *
 * which is EXACT for M/M/s/r+M (eq. 7.12), i.e. for `qsys_erlanga`.
 *
 * ARITHMETIC. The birth-death leg and every moment are field operations, so the
 * hazard and exponential forms instantiate at T = Rational. The ccdf form needs
 * a logarithm and the waiting-time cdfs need a numerical Laplace inversion, so
 * the first throws and the second is skipped unless T is double.
 *
 * DIVERGENCE from the printed eqs. (3.5)-(3.6): they read
 * delta_j = int_{(j-1)/lambda}^{j/lambda} h(t) dt and Delta_k = -log F^c(k/lambda),
 * which are cumulative hazards, i.e. dimensionless, while delta and Delta are
 * rates everywhere else in the paper. They are the AVERAGE hazard over an
 * interval of length 1/lambda, so the factor lambda is missing. Restoring it
 * makes the ccdf form reduce to the exact Erlang A rates under exponential
 * patience, which the paper states this approximation does (eq. 7.12).
 *
 * Reference: W. Whitt (2005). Engineering solution of a basic call-center model.
 * Management Science 51(2), 221-235.
 */

#include <algorithm>
#include <complex>
#include <cstddef>
#include <functional>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/lti/laplace_invert.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Which of the three accepted descriptions of the patience law is carried. */
enum class PatienceForm { Exponential, Hazard, Ccdf };

/**
 * The patience (time-to-abandon) law, in the three forms the algorithm accepts.
 *
 * The engineering solution needs it only through its hazard rate near the
 * origin, so these three carry exactly that.
 */
template <class T>
struct Patience {
    PatienceForm form = PatienceForm::Exponential;
    T theta = num_traits<T>::from_int(0);
    std::function<T(const T&)> fun;

    /** Exponential patience of rate theta, h(t) = theta, the Erlang A case. */
    static Patience exponential(const T& theta) {
        if (theta < num_traits<T>::from_int(0))
            throw InputError("Patience::exponential: the rate theta must be non-negative");
        Patience p;
        p.form = PatienceForm::Exponential;
        p.theta = theta;
        return p;
    }

    /** Patience given by its hazard rate h = f/(1-F). */
    static Patience hazard(std::function<T(const T&)> h) {
        Patience p;
        p.form = PatienceForm::Hazard;
        p.fun = std::move(h);
        return p;
    }

    /** Patience given by its complementary cdf G(t) = 1-F(t). */
    static Patience ccdf(std::function<T(const T&)> g) {
        Patience p;
        p.form = PatienceForm::Ccdf;
        p.fun = std::move(g);
        return p;
    }
};

/** Steady-state measures of a multiserver queue with customer abandonment. */
template <class T>
struct QsysAbandonResult {
    std::vector<T> queueLengthDist;   ///< P(N = k), k = 0..s+r
    T probLoss;                       ///< P(an arrival is blocked), 0 when r is infinite
    T probNoWait;                     ///< P(W = 0) among entering customers
    T probServed;                     ///< P(S), eventually served
    T probAbandon;                    ///< P(A) = 1 - P(S)
    T meanNumber;                     ///< E[N]
    T varNumber;                      ///< Var[N]
    T meanQueueLength;                ///< E[Q], Q = (N-s)^+
    T varQueueLength;                 ///< Var[Q]
    T utilization;                    ///< E[min(N,s)]/s
    T throughput;                     ///< lambda(1-P_loss)P(S)
    T abandonRate;                    ///< lambda(1-P_loss)P(A)
    T meanWaitServed;                 ///< E[W|S]
    T varWaitServed;                  ///< Var[W|S]
    T meanWaitAbandon;                ///< E[W|A]
    T varWaitAbandon;                 ///< Var[W|A]
    T meanWait;                       ///< E[W] over entering customers
    T secondMomentWait;               ///< E[W^2] over entering customers
    std::vector<T> abandonRates;      ///< delta_j
    std::vector<T> totalAbandonRates; ///< Delta_k, index k, Delta_0 = 0
    std::size_t numWaitingSpaces = 0; ///< waiting spaces used, the truncation when r is infinite
    bool exponentialPatience = false; ///< whether the answer is exact
    T patienceRate;                   ///< the exponential rate, unset otherwise
    std::vector<double> waitPoints;        ///< times of the cdfs, empty when not requested
    std::vector<double> cdfWaitServed;     ///< P(W <= t | S)
    std::vector<double> cdfWaitAbandon;    ///< P(W <= t | A)
    std::vector<double> cdfWait;           ///< P(W <= t)
};

/** Options of `qsys_mgisrgi_whitt`, all with the MATLAB defaults. */
struct MgisrgiOptions {
    std::vector<double> wPoints;          ///< times for the waiting-time cdfs
    std::size_t maxQueue = 100000;        ///< truncation level used when r is infinite
    double tol = 1e-14;                   ///< relative tail tolerance for that truncation
    std::string invMethod = "euler";      ///< Laplace inversion method
    std::size_t invN = 41;                ///< number of inversion nodes
};

namespace detail {

/** One step of eqs. (3.3)-(3.4) (hazard form) or (3.5)-(3.6) (ccdf form). */
template <class T>
void mgisrgi_rate_step(std::size_t j, const T& lambda, const T& delta_prev,
                       const Patience<T>& patience, T& delta_j, T& delta_tot) {
    const T t = num_traits<T>::from_int(static_cast<long>(j)) / lambda;
    if (patience.form == PatienceForm::Ccdf) {
        if constexpr (num_traits<T>::has_transcendental) {
            const T g = patience.fun(t);
            if (g <= num_traits<T>::from_int(0))
                throw InputError(
                    "qsys_mgisrgi_whitt: the patience ccdf vanishes, so every customer has "
                    "abandoned by then; supply a hazard instead");
            using std::log;
            delta_tot = -lambda * log(g);
            delta_j = delta_tot - delta_prev;
        } else {
            throw InputError(
                "qsys_mgisrgi_whitt: the ccdf form of the patience law needs a logarithm, "
                "which the exact arithmetic does not have; supply a hazard instead");
        }
    } else {
        delta_j = patience.form == PatienceForm::Exponential ? patience.theta : patience.fun(t);
        delta_tot = delta_prev + delta_j;
    }
    if (delta_j < num_traits<T>::from_int(0))
        throw InputError("qsys_mgisrgi_whitt: the patience law produced a negative abandonment rate");
}

/**
 * Eqs. (7.10)-(7.11): with k waiting, the total departure rate before the jth
 * departure epoch is s mu + Delta_k - Delta_{j-1}, of which delta_j is the share
 * belonging to the customer of interest. Fills rate, phi and the survival
 * product prod_{l<j}(1-phi_k(l)).
 */
template <class T>
void mgisrgi_kernel(std::size_t k, const T& smu, const std::vector<T>& dlt,
                    const std::vector<T>& delta, std::vector<T>& rate, std::vector<T>& phi,
                    std::vector<T>& surv) {
    const T one = num_traits<T>::from_int(1);
    rate.resize(k);
    phi.resize(k);
    surv.resize(k);
    T running = one;
    for (std::size_t j = 1; j <= k; ++j) {
        rate[j - 1] = smu + dlt[k] - dlt[j - 1];
        phi[j - 1] = delta[j - 1] / rate[j - 1];
        surv[j - 1] = running;
        running *= (one - phi[j - 1]);
    }
}

/** A conditional moment is 0/0 when the conditioning event cannot happen. */
template <class T>
T mgisrgi_ratio(const T& num, const T& den) {
    return den <= num_traits<T>::from_int(0) ? num_traits<T>::from_int(0) : T(num / den);
}

}  // namespace detail

/**
 * @param lambda   arrival rate
 * @param mu       service rate of one server, the reciprocal of the mean service time
 * @param s        number of servers, s >= 1
 * @param r        extra waiting spaces; infinity for an unbounded queue
 * @param patience the patience law
 * @param opts     cdf times, truncation controls and inversion settings
 */
template <class T>
QsysAbandonResult<T> qsys_mgisrgi_whitt(const T& lambda, const T& mu, unsigned s, double r,
                                        const Patience<T>& patience,
                                        const MgisrgiOptions& opts = MgisrgiOptions()) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (lambda <= zero) throw InputError("qsys_mgisrgi_whitt: the arrival rate lambda must be positive");
    if (mu <= zero) throw InputError("qsys_mgisrgi_whitt: the service rate mu must be positive");
    if (s < 1) throw InputError("qsys_mgisrgi_whitt: the number of servers s must be at least 1");
    if (r < 0)
        throw InputError("qsys_mgisrgi_whitt: the number of extra waiting spaces r must be non-negative");

    const bool finiteR = std::isfinite(r);
    std::size_t rr = finiteR ? static_cast<std::size_t>(r + 0.5) : opts.maxQueue;
    const T smu = num_traits<T>::from_int(static_cast<long>(s)) * mu;

    // The birth-death recursion of eqs. (7.4)-(7.7), unnormalized with x_s = 1.
    std::vector<T> xUp(rr + 1, zero), dlt(rr + 1, zero), delta(rr, zero);
    xUp[0] = one;
    std::size_t kUsed = rr;
    T peak = one;
    const T tolT = num_traits<T>::from_double(opts.tol);
    for (std::size_t k = 0; k < rr; ++k) {
        const std::size_t j = k + 1;
        detail::mgisrgi_rate_step(j, lambda, dlt[j - 1], patience, delta[j - 1], dlt[j]);
        xUp[k + 1] = lambda * xUp[k] / (smu + dlt[j]);
        if (xUp[k + 1] > peak) peak = xUp[k + 1];
        if (!finiteR && xUp[k + 1] < tolT * peak && k >= 1) {
            kUsed = j;
            break;
        }
    }
    if (!finiteR) {
        if (kUsed == rr && rr > 0)
            throw InputError("qsys_mgisrgi_whitt: the queue-length tail is still sizeable at the "
                             "truncation level; with r = Inf the patience law must make the chain "
                             "ergodic (raise maxQueue if the model is genuinely that large)");
        xUp.resize(kUsed + 1);
        dlt.resize(kUsed + 1);
        delta.resize(kUsed);
        rr = kUsed;
    }

    // The downward leg, eq. (7.5), over the states where not all servers are busy.
    std::vector<T> x(s + rr + 1, zero);
    T xk = one;
    for (unsigned k = s; k >= 1; --k) {
        xk = num_traits<T>::from_int(static_cast<long>(k)) * mu * xk / lambda;
        x[k - 1] = xk;
    }
    for (std::size_t k = 0; k <= rr; ++k) x[s + k] = xUp[k];

    T total = zero;
    for (const T& v : x) total += v;
    QsysAbandonResult<T> res;
    res.queueLengthDist.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) res.queueLengthDist[i] = x[i] / total;
    res.probLoss = finiteR ? res.queueLengthDist.back() : zero;
    std::vector<T> pa(res.queueLengthDist.size());
    for (std::size_t i = 0; i < pa.size(); ++i)
        pa[i] = res.queueLengthDist[i] / (one - res.probLoss);   // eq. (7.8), seen by an ENTERING customer

    T meanNumber = zero, meanQueue = zero, busy = zero;
    for (std::size_t k = 0; k < res.queueLengthDist.size(); ++k) {
        const T pk = res.queueLengthDist[k];
        const T kT = num_traits<T>::from_int(static_cast<long>(k));
        meanNumber += kT * pk;
        if (k > s) meanQueue += num_traits<T>::from_int(static_cast<long>(k - s)) * pk;
        busy += num_traits<T>::from_int(static_cast<long>(std::min<std::size_t>(k, s))) * pk;
    }
    T varNumber = zero, varQueue = zero;
    for (std::size_t k = 0; k < res.queueLengthDist.size(); ++k) {
        const T pk = res.queueLengthDist[k];
        const T dN = num_traits<T>::from_int(static_cast<long>(k)) - meanNumber;
        const T q = k > s ? num_traits<T>::from_int(static_cast<long>(k - s)) : zero;
        varNumber += dN * dN * pk;
        varQueue += (q - meanQueue) * (q - meanQueue) * pk;
    }

    T probNoWait = zero;                       // eq. (7.9), states 0..s-1
    for (unsigned k = 0; k < s; ++k) probNoWait += pa[k];

    std::vector<T> sigma(rr, zero), mSum(rr, zero), vSum(rr, zero), ewa1(rr, zero), ewa2(rr, zero);
    std::vector<std::vector<T>> rateK(rr), phiK(rr), survK(rr);
    for (std::size_t k = 1; k <= rr; ++k) {
        detail::mgisrgi_kernel(k, smu, dlt, delta, rateK[k - 1], phiK[k - 1], survK[k - 1]);
        T prod = one, sm = zero, sv = zero, cumM = zero, cumV = zero, e1 = zero, e2 = zero;
        for (std::size_t j = 0; j < k; ++j) {
            const T m = one / rateK[k - 1][j];
            prod *= (one - phiK[k - 1][j]);
            sm += m;
            sv += m * m;
            // Eqs. (7.28)-(7.29): abandoning at the jth departure epoch costs the
            // sum of the first j interdeparture times.
            cumM += m;
            cumV += m * m;
            const T w = survK[k - 1][j] * phiK[k - 1][j];
            e1 += w * cumM;
            e2 += w * (cumV + cumM * cumM);
        }
        sigma[k - 1] = prod;
        mSum[k - 1] = sm;
        vSum[k - 1] = sv;
        ewa1[k - 1] = e1;
        ewa2[k - 1] = e2;
    }

    // Finding s+k in system puts the arrival in position k+1, so the weights are
    // pa_{s+k} for k = 0..r-1.
    std::vector<T> wArr(rr, zero);
    for (std::size_t k = 0; k < rr; ++k) wArr[k] = pa[s + k];
    T probServed = probNoWait, ews1 = zero, ews2 = zero, ewa1Tot = zero, ewa2Tot = zero;
    for (std::size_t k = 0; k < rr; ++k) {
        probServed += wArr[k] * sigma[k];
        ews1 += wArr[k] * sigma[k] * mSum[k];                            // eq. (7.16)
        ews2 += wArr[k] * sigma[k] * (vSum[k] + mSum[k] * mSum[k]);      // eq. (7.17)
        ewa1Tot += wArr[k] * ewa1[k];                                    // eq. (7.26)
        ewa2Tot += wArr[k] * ewa2[k];                                    // eq. (7.27)
    }

    res.probNoWait = probNoWait;
    res.probServed = probServed;
    res.probAbandon = one - probServed;
    res.meanNumber = meanNumber;
    res.varNumber = varNumber;
    res.meanQueueLength = meanQueue;
    res.varQueueLength = varQueue;
    res.utilization = busy / num_traits<T>::from_int(static_cast<long>(s));
    res.throughput = lambda * (one - res.probLoss) * probServed;
    res.abandonRate = lambda * (one - res.probLoss) * res.probAbandon;
    res.meanWaitServed = detail::mgisrgi_ratio(ews1, probServed);
    res.varWaitServed = detail::mgisrgi_ratio(ews2, probServed) -
                        res.meanWaitServed * res.meanWaitServed;
    if (res.varWaitServed < zero) res.varWaitServed = zero;
    res.meanWaitAbandon = detail::mgisrgi_ratio(ewa1Tot, res.probAbandon);
    res.varWaitAbandon = detail::mgisrgi_ratio(ewa2Tot, res.probAbandon) -
                         res.meanWaitAbandon * res.meanWaitAbandon;
    if (res.varWaitAbandon < zero) res.varWaitAbandon = zero;
    res.meanWait = ews1 + ewa1Tot;
    res.secondMomentWait = ews2 + ewa2Tot;
    res.abandonRates = delta;
    res.totalAbandonRates = dlt;
    res.numWaitingSpaces = rr;
    res.exponentialPatience = patience.form == PatienceForm::Exponential;
    res.patienceRate = patience.theta;

    if (!opts.wPoints.empty()) {
        if constexpr (std::is_same_v<T, double>) {
            // Eqs. (7.22)-(7.23) served, (7.32)-(7.33) abandoning. Both fold the
            // same kernel: the wait is a sum of exponentials with rates 1/m_k(j),
            // truncated at the departure epoch that serves or loses the customer.
            auto transform = [&](const lti::Cplx& z, bool served) {
                lti::Cplx val(0.0, 0.0);
                for (std::size_t k = 1; k <= rr; ++k) {
                    lti::Cplx chain(1.0, 0.0);
                    for (std::size_t j = 0; j < k; ++j) {
                        chain *= rateK[k - 1][j] / (rateK[k - 1][j] + z);
                        if (!served) val += chain * (wArr[k - 1] * survK[k - 1][j] * phiK[k - 1][j]);
                    }
                    if (served) val += chain * (wArr[k - 1] * sigma[k - 1]);
                }
                return val;
            };
            const lti::LaplaceMethod method = lti::laplace_method(opts.invMethod);
            const double capS = std::max(probServed - probNoWait, 0.0);
            const double capA = std::max(res.probAbandon, 0.0);
            res.waitPoints = opts.wPoints;
            res.cdfWaitServed.resize(opts.wPoints.size());
            res.cdfWaitAbandon.resize(opts.wPoints.size());
            res.cdfWait.resize(opts.wPoints.size());
            for (std::size_t i = 0; i < opts.wPoints.size(); ++i) {
                double fs = lti::laplace_invert(
                    [&](const lti::Cplx& z) { return transform(z, true) / z; }, opts.wPoints[i],
                    method, opts.invN);
                double fa = lti::laplace_invert(
                    [&](const lti::Cplx& z) { return transform(z, false) / z; }, opts.wPoints[i],
                    method, opts.invN);
                fs = std::min(std::max(fs, 0.0), capS);
                fa = std::min(std::max(fa, 0.0), capA);
                res.cdfWaitServed[i] = (probNoWait + fs) / std::max(probServed, 1e-300);
                res.cdfWaitAbandon[i] = fa / std::max(res.probAbandon, 1e-300);
                res.cdfWait[i] = probNoWait + fs + fa;
            }
        } else {
            throw InputError("qsys_mgisrgi_whitt: the waiting-time cdfs need a numerical Laplace "
                             "inversion, which only the double instantiation carries");
        }
    }
    return res;
}

/**
 * Exact analysis of the Erlang A model M/M/s/r+M.
 *
 * The number in system is the birth-and-death process with death rate
 * min(k,s) mu + (k-s)^+ theta, so every measure is exact: this is the case in
 * which the approximation above reproduces the model (eq. 7.12). theta = 0
 * recovers M/M/s/r, and then a finite r is required whenever lambda >= s mu.
 *
 * Port of matlab/src/api/qsys/qsys_erlanga.m.
 *
 * @param lambda arrival rate
 * @param mu     service rate of one server
 * @param theta  abandonment rate of a waiting customer
 * @param s      number of servers
 * @param r      extra waiting spaces; infinity for an unbounded queue
 * @param opts   cdf times, truncation controls and inversion settings
 */
template <class T>
QsysAbandonResult<T> qsys_erlanga(const T& lambda, const T& mu, const T& theta, unsigned s,
                                  double r = std::numeric_limits<double>::infinity(),
                                  const MgisrgiOptions& opts = MgisrgiOptions()) {
    if (theta <= num_traits<T>::from_int(0) && !std::isfinite(r) &&
        lambda >= num_traits<T>::from_int(static_cast<long>(s)) * mu)
        throw InputError("qsys_erlanga: without abandonment (theta = 0) and with an infinite "
                         "waiting room the queue is unstable at lambda >= s*mu; give a finite r "
                         "or a positive theta");
    return qsys_mgisrgi_whitt(lambda, mu, s, r, Patience<T>::exponential(theta), opts);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MGISRGI_WHITT_H
