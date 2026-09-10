/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MCMC_H
#define LINE_API_PFQN_MCMC_H

/**
 * Chen-O'Cinneide REGULARIZATION: a Markov chain Monte Carlo estimator of the
 * class throughputs X(r) = G(N-e_r)/G(N) and of the mean queue lengths Q(i,r) of
 * a CLOSED multiclass product-form (BCMP, no type changes) network.
 *
 *   W. Chen, C. A. O'Cinneide, "Towards a Polynomial-Time Randomized Algorithm
 *   for Closed Product-Form Networks", ACM TOMACS 8(3):227-253, 1998.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mcmc.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/nc/Pfqn_mcmc.java.
 *
 * The three steps of the paper are:
 *
 *  I.  CONSTRUCT THE REGULARIZED NETWORK. Write rho(i,r) for the surrogate
 *      traffic intensity of class r at station i -- here the service demand,
 *      since rho = lambda/mu is a visit ratio over a service rate -- and
 *      rho(r) = sum_i rho(i,r). The regularized network has the same stations,
 *      classes and populations, UNIT service rates at every station, the
 *      processor-sharing discipline, and a routing matrix that depends on the
 *      destination only,
 *
 *          P*(i->m | class r) = rho(m,r)/rho(r).
 *
 *      By Theorem 2.1 it is a REVERSIBLE chain with the SAME steady-state
 *      distribution as the original network, and its throughputs satisfy
 *      Theta*(r) = rho(r)*Theta(r).
 *
 *  II. SIMULATE IT at service-completion epochs. With Y(i,r) the number of
 *      class-r jobs at station i, Y(i) their total and Psi_i(k) = min(s_i,k) the
 *      number of busy servers,
 *
 *          r(i,r) = Y(i,r)/Y(i) * Psi_i(Y(i)),   r(r) = sum_i r(i,r),
 *          r      = sum_i Psi_i(Y(i)),
 *
 *      the next completion is of class r at station i with probability
 *      r(i,r)/r, and the conditional expected time to it is 1/r. Equation (10)
 *      of the paper is the holding-time weighted ratio estimator
 *
 *          Theta*(r) = sum_t r(r,t)/r(t)  /  sum_t 1/r(t),
 *
 *      and the same weights give the time-average queue lengths, which need no
 *      transformation at all because the two networks share their steady state.
 *
 * III. TRANSFORM BACK: X(r) = Theta*(r)/rho(r).
 *
 * Because P* forgets the station of origin and every station serves at unit
 * rate, the regularized chain has neither the slowly mixing routing chain nor
 * the customer-trapping slow station that make the original chain converge
 * slowly. The paper proves O(N^2*M^3) mixing in two special cases (Section 4)
 * and reports the general behaviour experimentally (Section 5).
 *
 * Delay (infinite-server) demand enters as ONE extra station with s = infinity
 * and demand Z. Aggregating infinite-server stations that way is exact in the
 * product form, since their joint term is multinomial in the per-class totals.
 *
 * THERE IS NO NORMALIZING CONSTANT HERE, and that is the single most important
 * thing to know about the method: the estimator is a ratio, G itself never
 * appears. `pfqn_nc`'s Mcmc arm supplies an BLE lG alongside, which is not
 * part of the paper and cancels out of every mean value.
 *
 * Confidence: the run is split into non-overlapping batches (Schmeiser 1982, 30
 * by default, the count used in the tables of the paper), the batch means of the
 * ratio estimator give a standard error, and the intervals are the paper's
 * two-sigma ones. The estimator is a ratio of correlated averages, so it carries
 * an O(1/samples) bias on top of the initialization bias; the paper ignores
 * both, this port additionally discards a warm-up fraction (10% by default).
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The value is a random variable and the
 * batch-means interval needs a square root, so the estimator is meaningless --
 * not merely inaccurate -- in an exact field, exactly as the rest of the
 * `pfqn_nc` estimator ladder is.
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Sample standard deviation over the batch means, normalized by n-1 as MATLAB's std is. */
inline double mcmc_sample_std(const std::vector<double>& v) {
    const std::size_t n = v.size();
    if (n < 2) return 0.0;
    double mean = 0.0;
    for (std::size_t i = 0; i < n; ++i) mean += v[i];
    mean /= static_cast<double>(n);
    double acc = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = v[i] - mean;
        acc += d * d;
    }
    return std::sqrt(acc / static_cast<double>(n - 1));
}

}  // namespace detail

/** Estimates of pfqn_mcmc together with their batch-means intervals. */
template <class T>
struct McmcResult {
    std::vector<T> X;   ///< (R) throughput estimates G(N-e_r)/G(N)
    Matrix<T> Q;        ///< (M x R) mean queue lengths at the queueing stations
    std::vector<T> Xse; ///< (R) batch-means standard error of X
    std::vector<T> Xlo; ///< (R) lower end of the two-sigma interval for X
    std::vector<T> Xhi; ///< (R) upper end of the two-sigma interval for X
    Matrix<T> Qse;      ///< (M x R) batch-means standard error of Q
    Matrix<T> Qlo;      ///< (M x R) lower end of the two-sigma interval for Q
    Matrix<T> Qhi;      ///< (M x R) upper end of the two-sigma interval for Q
    std::size_t batches = 0;  ///< batches the run was split into
    std::size_t samples = 0;  ///< completions simulated after warm-up
    std::size_t burnin = 0;   ///< completions discarded as warm-up
};

/** Schmeiser (1982), the batch count used in the tables of the paper. */
inline constexpr std::size_t MCMC_DEFAULT_BATCHES = 30;
/** Warm-up fraction discarded before accumulation starts. */
inline constexpr double MCMC_DEFAULT_BURNIN = 0.1;

/**
 * @param L        (M x R) per-class service demands at the M queueing stations
 * @param N        (R) closed population vector; finite and integer
 * @param Z        (R) aggregated think times; empty for a model with no delay
 * @param s        (M) servers per station, infinite for an infinite server;
 *                 empty means all stations single-server
 * @param samples  service completions to simulate after warm-up
 * @param nbatches batches the run is split into for the confidence intervals
 * @param burnin   warm-up fraction discarded before accumulation starts
 * @param rng      explicit generator, advanced by the call
 */
template <class T>
McmcResult<T> pfqn_mcmc(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                        const std::vector<double>& s, std::size_t samples, std::size_t nbatches,
                        double burnin, McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_mcmc requires transcendental arithmetic: it is a Monte Carlo estimator "
                  "whose value is a random variable and whose batch-means interval needs a "
                  "square root");

    const std::size_t M = L.empty() ? 0 : L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_mcmc: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_mcmc: Z has the wrong length");
    if (!s.empty() && s.size() != M)
        throw InputError("pfqn_mcmc: the server count vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    McmcResult<T> res;
    res.X.assign(R, zero);
    res.Xse.assign(R, zero);
    res.Xlo.assign(R, zero);
    res.Xhi.assign(R, zero);
    res.Q = Matrix<T>(M, R, zero);
    res.Qse = Matrix<T>(M, R, zero);
    res.Qlo = Matrix<T>(M, R, zero);
    res.Qhi = Matrix<T>(M, R, zero);

    long Ntot = 0;
    for (std::size_t r = 0; r < R; ++r) {
        // The chain lives on the integer lattice sum_i Y(i,r) = N(r); the caller has
        // already resolved the population to an int, so a fractional one cannot reach
        // here, but a NEGATIVE entry is the open-class marker and has no state space.
        if (N[r] < 0)
            throw InputError("pfqn_mcmc requires a closed model, but the population vector "
                             "marks an open class");
        Ntot += N[r];
    }
    if (Ntot == 0) return res;

    // ---- Step I: the regularized network ---------------------------------------------
    // Only the surrogate traffic intensities rho(i,r) enter the product form, and scaling
    // a whole class column by a constant leaves the steady-state distribution unchanged,
    // so the demands are used as they are.
    bool hasDelay = false;
    for (std::size_t r = 0; r < R && !Z.empty(); ++r)
        if (num_traits<T>::to_double(Z[r]) > 0.0) hasDelay = true;
    const std::size_t Mx = hasDelay ? M + 1 : M;

    Matrix<T> rho(Mx, R, zero);
    std::vector<double> svec(Mx, 1.0);
    for (std::size_t i = 0; i < M; ++i) {
        svec[i] = s.empty() ? 1.0 : s[i];
        for (std::size_t r = 0; r < R; ++r) {
            const double v = num_traits<T>::to_double(L(i, r));
            rho(i, r) = (std::isfinite(v) && v > 0.0) ? L(i, r) : zero;
        }
    }
    if (hasDelay) {
        svec[M] = std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < R; ++r) rho(M, r) = Z[r];
    }

    std::vector<T> rhoTot(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        for (std::size_t i = 0; i < Mx; ++i) rhoTot[r] += rho(i, r);
        if (N[r] > 0 && num_traits<T>::to_double(rhoTot[r]) <= 0.0)
            throw InputError("pfqn_mcmc: a class has a positive population but no demand "
                             "anywhere in the network");
    }

    // Routing of the regularized network, P*(m|r) = rho(m,r)/rho(r), held as one column of
    // cumulative probabilities per class.
    Matrix<double> Pstar(Mx, R, 0.0);
    Matrix<double> cumP(Mx, R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        const double den = std::max(num_traits<T>::to_double(rhoTot[r]),
                                    std::numeric_limits<double>::min());
        double acc = 0.0;
        for (std::size_t i = 0; i < Mx; ++i) {
            Pstar(i, r) = num_traits<T>::to_double(rho(i, r)) / den;
            acc += Pstar(i, r);
            cumP(i, r) = acc;
        }
        cumP(Mx - 1, r) = 1.0;  // guard the last bin against a floating-point shortfall
    }

    if (nbatches == 0) nbatches = 1;
    if (samples == 0) samples = 1;
    if (!(burnin >= 0.0)) burnin = 0.0;
    if (burnin > 0.9) burnin = 0.9;
    const std::size_t batchLen = std::max<std::size_t>(1, samples / nbatches);
    samples = batchLen * nbatches;
    const std::size_t nburn = static_cast<std::size_t>(std::llround(burnin * samples));

    // Initial state: spread each class over the stations it can occupy in the proportions
    // P*(.|r), by largest remainder. That is the marginal the regularized network would
    // have with no queueing, so it costs nothing and starts the chain far closer to
    // stationarity than a single-station state.
    std::vector<std::vector<long> > Y(Mx, std::vector<long>(R, 0));
    std::vector<double> Ytot(Mx, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] == 0) continue;
        std::vector<double> target(Mx, 0.0);
        long placed = 0;
        for (std::size_t i = 0; i < Mx; ++i) {
            target[i] = N[r] * Pstar(i, r);
            Y[i][r] = static_cast<long>(std::floor(target[i]));
            placed += Y[i][r];
        }
        // Largest remainder: hand each remaining unit to the biggest residue still
        // outstanding. Once a station is topped up its residue target-Y drops by one, so
        // it cannot win twice, which is the 'descend' sort of the reference.
        for (long k = 0; k < N[r] - placed; ++k) {
            std::size_t best = 0;
            double bestRem = -std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < Mx; ++i) {
                const double rem = target[i] - static_cast<double>(Y[i][r]);
                if (rem > bestRem) {
                    bestRem = rem;
                    best = i;
                }
            }
            ++Y[best][r];
        }
    }
    for (std::size_t i = 0; i < Mx; ++i) {
        double tot = 0.0;
        for (std::size_t r = 0; r < R; ++r) tot += static_cast<double>(Y[i][r]);
        Ytot[i] = tot;
    }

    // ---- Step II: simulate at service-completion epochs -------------------------------
    Matrix<double> xnum(nbatches, R, 0.0);           // sum_t r(r,t)/r(t) within the batch
    std::vector<Matrix<double> > qnum(nbatches, Matrix<double>(Mx, R, 0.0));
    std::vector<double> den(nbatches, 0.0);          // sum_t 1/r(t) within the batch
    std::vector<double> Psi(Mx, 0.0), cPsi(Mx, 0.0), rvec(R, 0.0), cY(R, 0.0);
    double w = 0.0;
    bool stale = true;
    const std::size_t horizon = nburn + samples;
    for (std::size_t t = 1; t <= horizon; ++t) {
        if (stale) {
            // (8)-(9): busy servers, per-class completion rates, total rate
            double totPsi = 0.0;
            for (std::size_t i = 0; i < Mx; ++i) {
                Psi[i] = std::min(svec[i], Ytot[i]);
                totPsi += Psi[i];
                cPsi[i] = totPsi;
            }
            for (std::size_t r = 0; r < R; ++r) rvec[r] = 0.0;
            for (std::size_t i = 0; i < Mx; ++i) {
                if (Ytot[i] <= 0.0) continue;
                const double rw = Psi[i] / Ytot[i];
                for (std::size_t r = 0; r < R; ++r)
                    if (Y[i][r] != 0) rvec[r] += rw * static_cast<double>(Y[i][r]);
            }
            w = 1.0 / totPsi;
            stale = false;
        }
        if (t > nburn) {
            const std::size_t b = (t - nburn - 1) / batchLen;
            den[b] += w;
            for (std::size_t r = 0; r < R; ++r) xnum(b, r) += w * rvec[r];
            Matrix<double>& qb = qnum[b];
            for (std::size_t i = 0; i < Mx; ++i)
                for (std::size_t r = 0; r < R; ++r)
                    if (Y[i][r] != 0) qb(i, r) += w * static_cast<double>(Y[i][r]);
        }
        // Pick the completing station with probability Psi(i)/r, then the completing class
        // within it with probability Y(i,r)/Y(i); the product is the r(i,r)/r of the
        // paper, since sum_r Y(i,r)/Y(i)*Psi(i) = Psi(i).
        std::size_t i = Mx;
        {
            const double u = mc_uniform01(rng) * cPsi[Mx - 1];
            for (std::size_t k = 0; k < Mx; ++k)
                if (cPsi[k] >= u) {
                    i = k;
                    break;
                }
            if (i == Mx)
                for (std::size_t k = Mx; k-- > 0;)
                    if (Psi[k] > 0.0) {
                        i = k;
                        break;
                    }
        }
        double acc = 0.0;
        for (std::size_t r = 0; r < R; ++r) {
            acc += static_cast<double>(Y[i][r]);
            cY[r] = acc;
        }
        std::size_t cls = R;
        {
            const double u = mc_uniform01(rng) * cY[R - 1];
            for (std::size_t k = 0; k < R; ++k)
                if (cY[k] >= u) {
                    cls = k;
                    break;
                }
            if (cls == R)
                for (std::size_t k = R; k-- > 0;)
                    if (Y[i][k] > 0) {
                        cls = k;
                        break;
                    }
        }
        // Route it. A self-transition leaves the state, hence the rates and the weight,
        // unchanged: skipping the recomputation is the saving described at the end of
        // Section 2 of the paper.
        const double u = mc_uniform01(rng);
        std::size_t m = Mx;
        for (std::size_t k = 0; k < Mx; ++k)
            if (cumP(k, cls) >= u) {
                m = k;
                break;
            }
        if (m < Mx && m != i) {
            --Y[i][cls];
            ++Y[m][cls];
            Ytot[i] -= 1.0;
            Ytot[m] += 1.0;
            stale = true;
        }
    }

    // ---- Step III: back to the original network ---------------------------------------
    // Theta(r) = Theta*(r)/rho(r) by (7) and (11); the queue lengths transfer unchanged,
    // the two networks sharing their steady-state distribution.
    double denTot = 0.0;
    for (std::size_t b = 0; b < nbatches; ++b) denTot += den[b];
    for (std::size_t r = 0; r < R; ++r) {
        double num = 0.0;
        for (std::size_t b = 0; b < nbatches; ++b) num += xnum(b, r);
        res.X[r] = num_traits<T>::from_double((num / denTot) /
                                              num_traits<T>::to_double(rhoTot[r]));
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            double num = 0.0;
            for (std::size_t b = 0; b < nbatches; ++b) num += qnum[b](i, r);
            res.Q(i, r) = num_traits<T>::from_double(num / denTot);
        }

    // Batch-means standard error and the two-sigma interval of the paper.
    if (nbatches > 1) {
        std::vector<double> v(nbatches, 0.0);
        for (std::size_t r = 0; r < R; ++r) {
            const double rt = num_traits<T>::to_double(rhoTot[r]);
            for (std::size_t b = 0; b < nbatches; ++b) v[b] = (xnum(b, r) / den[b]) / rt;
            res.Xse[r] = num_traits<T>::from_double(detail::mcmc_sample_std(v) /
                                                    std::sqrt(static_cast<double>(nbatches)));
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t b = 0; b < nbatches; ++b) v[b] = qnum[b](i, r) / den[b];
                res.Qse(i, r) = num_traits<T>::from_double(
                    detail::mcmc_sample_std(v) / std::sqrt(static_cast<double>(nbatches)));
            }
    }
    const T two = num_traits<T>::from_int(2);
    for (std::size_t r = 0; r < R; ++r) {
        res.Xlo[r] = T(res.X[r] - two * res.Xse[r]);
        res.Xhi[r] = T(res.X[r] + two * res.Xse[r]);
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            res.Qlo(i, r) = T(res.Q(i, r) - two * res.Qse(i, r));
            res.Qhi(i, r) = T(res.Q(i, r) + two * res.Qse(i, r));
        }
    res.batches = nbatches;
    res.samples = samples;
    res.burnin = nburn;
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MCMC_H
