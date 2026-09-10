/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_GIBBS_H
#define LINE_API_INFER_INFER_GIBBS_H

/**
 * Gibbs sampling demand estimator for a closed delay-plus-queue model.
 *
 * Templated port of matlab/src/api/infer/infer_gibbs.m. No JAR counterpart.
 *
 * The observed data are per-class arrival and departure traces at one station.
 * The estimator replays them into a continuous-time record of the population
 * vector, turns that record into an empirical distribution over states,
 * samples a test set from it, and then draws the per-class demands from their
 * conditional posteriors one coordinate at a time. The posterior of a single
 * demand is a slice: on a grid of candidate values,
 *
 *   log p(theta_j) = (sum over the test set of the class-j queue count)
 *                    log theta_j - n_test log G(theta),
 *
 * i.e. the product-form likelihood with the normalizing constant carried
 * along. G is never evaluated directly. It is propagated along the grid with
 * the exact derivative identity of the closed product-form constant,
 * d log G / d theta_j = Q_j / theta_j, integrated one grid step at a time with
 * Q from Bard-Schweitzer AMVA (the 'TE' method of the MATLAB). That is why
 * every grid step costs one pfqn_bs solve, warm started from the previous
 * point, and why the walk starts at the current theta with the running log G
 * the caller carries between coordinate updates.
 *
 * WHAT IS AND IS NOT PORTED. MATLAB's 'MCI' branch and its pdf_slice helper
 * are unreachable: alg is assigned 'TE' as a literal with no way in, and
 * pdf_slice is DECLARED with eleven parameters and CALLED with ten from the
 * one dead call site, so the MCI branch would raise on its first use. Dead,
 * broken code is not ported. Everything the 'TE' path executes is.
 *
 * The four sample budgets that MATLAB hard-codes (data_needed, the test set
 * size, the chain length, and the 50-sample convergence block) are options
 * here, defaulting to MATLAB's values. That is a strict superset: the default
 * construction reproduces the MATLAB exactly, and a caller that wants a short
 * chain no longer has to edit the source.
 *
 * RANDOMNESS. Two places draw: the test set is sampled from the empirical
 * state distribution by inversion, and each coordinate update is drawn from
 * its normalized slice. Both take an explicit McRng, so a run is reproducible
 * from its seed. MATLAB's Mersenne Twister stream is deliberately NOT
 * reproduced -- same algorithm, different stream -- so the two implementations
 * agree in distribution and on every deterministic intermediate, not sample by
 * sample.
 *
 * MATLAB's `eps` in the derivative denominators is the DOUBLE machine epsilon,
 * a fixed constant of the formula rather than a property of the working
 * arithmetic, so the port carries the literal 2^-52 into T instead of asking
 * the arithmetic for its own epsilon. Substituting Real50's epsilon would
 * change the first grid step, where theta = 0 and the constant is the entire
 * denominator.
 *
 * ARITHMETIC: logarithms throughout, plus a fixed-point AMVA iteration per
 * grid point, so the routine is gated on transcendental arithmetic and
 * registered for Double only. Real is NOT registered: log(0) at the first grid
 * point is a defined -infinity in IEEE double and the slice weight it produces
 * underflows to exactly zero, which is the behaviour the algorithm relies on;
 * that is a property of the double arithmetic MATLAB runs in, and the port
 * does not claim a high-precision instantiation it has not validated.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/** The two nodes of the model the estimator assumes: delay, then queue. */
static const std::size_t GIBBS_NNODES = 2;

/**
 * One class's trace. MATLAB reads these from rows 3, 4 and 6 of its 6 x (K+1)
 * cell array; the cell container is a MATLAB storage detail, so the port takes
 * the three sample vectors directly, as infer_get_qlen_arrival already does.
 */
template <class T>
struct GibbsTrace {
    std::vector<T> arrival_ms;  ///< data{3,k}, arrival times in MILLISECONDS
    std::vector<T> respt_s;     ///< data{4,k}, response times in SECONDS
    std::vector<T> think_obs;   ///< data{6,k}, the think-time normalizer
};

/** The empirical state distribution built from the traces. */
template <class T>
struct GibbsStateProbs {
    Matrix<long> states;    ///< (ns x K*GIBBS_NNODES) population vectors, node-major
    std::vector<T> prob;    ///< (ns) time fraction spent in each state
    std::vector<long> N;    ///< (K) per-class population
    std::vector<T> N0;      ///< (K) per-class mean number at the queue
};

/** MATLAB's hard-coded budgets, exposed with their MATLAB values as defaults. */
struct GibbsOptions {
    double tol = 1e-3;                     ///< grid step and convergence tolerance
    std::size_t data_needed = 200000;      ///< events kept from the end of the trace
    std::size_t likelihood_sample = 5000;  ///< test set size
    std::size_t nsamples = 2000;           ///< chain length
    std::size_t block = 50;                ///< samples per convergence block
};

/** The deterministic content of one coordinate update. */
template <class T>
struct GibbsSlice {
    std::vector<T> grid;   ///< candidate values, 0, interval, 2 interval, ...
    std::vector<T> logG;   ///< log normalizing constant along the grid
    std::vector<T> prob;   ///< normalized slice probability
    T range_size_dim;      ///< twice the value at which the slice mass closes
};

namespace detail {

/** MATLAB's eps, the double machine epsilon, as a constant of the formula. */
template <class T>
T gibbs_eps() {
    return num_traits<T>::from_double(2.220446049250313e-16);
}

/** Number of points of MATLAB's 0:step:limit, which counts floor(limit/step)+1. */
inline std::size_t colon_count(double limit, double step) {
    if (!(step > 0.0)) throw InputError("infer_gibbs: the grid step must be positive");
    if (limit < 0.0) return 0;
    return static_cast<std::size_t>(std::floor(limit / step + 1e-10)) + 1;
}

}  // namespace detail

/**
 * Empirical state distribution of the replayed traces (MATLAB's analyseData).
 *
 * Each class contributes one event per sample at the arrival time and one at
 * the arrival plus the response time, so the record is a birth-death walk of
 * the population vector between the delay node and the queue. The population
 * of each class is taken as the largest count the walk ever reaches, which is
 * what makes the delay-node counts non-negative once it is added back.
 *
 * The state is integer by construction, so it is carried as long and only the
 * holding times are in T: the aggregation over repeated states is then exact
 * in every arithmetic instead of accumulating a rounding per event.
 *
 * @param data        per-class traces
 * @param data_needed events kept from the END of the record, 0 for all
 */
template <class T>
GibbsStateProbs<T> gibbs_analyse_data(const std::vector<GibbsTrace<T>>& data,
                                      std::size_t data_needed) {
    const std::size_t K = data.size();
    if (K == 0) throw InputError("gibbs_analyse_data: no classes");

    // Event stream: an arrival and a departure per sample of every class.
    struct Ev {
        T t;
        std::size_t cls;
        int logger;  // 1 = into the queue, 2 = back to the delay
    };
    const T thousand = num_traits<T>::from_int(1000);
    std::vector<Ev> ev;
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t n = data[k].arrival_ms.size();
        if (data[k].respt_s.size() != n)
            throw InputError("gibbs_analyse_data: a class has mismatched sample counts");
        for (std::size_t i = 0; i < n; ++i) {
            Ev a;
            a.t = data[k].arrival_ms[i];
            a.cls = k;
            a.logger = 1;
            ev.push_back(a);
        }
        for (std::size_t i = 0; i < n; ++i) {
            Ev d;
            d.t = T(data[k].arrival_ms[i] + data[k].respt_s[i] * thousand);
            d.cls = k;
            d.logger = 2;
            ev.push_back(d);
        }
    }
    const std::size_t total = ev.size();
    if (total == 0) throw InputError("gibbs_analyse_data: empty traces");
    // MATLAB's sort is stable, and ties between an arrival and a departure at
    // the same instant change the state sequence, so the order matters.
    std::stable_sort(ev.begin(), ev.end(),
                     [](const Ev& a, const Ev& b) { return a.t < b.t; });

    // Replay. count(i) is the state that HOLDS from ev[i-1].t to ev[i].t.
    Matrix<long> count(total, K * GIBBS_NNODES, 0L);
    for (std::size_t i = 0; i + 1 < total; ++i) {
        for (std::size_t c = 0; c < K * GIBBS_NNODES; ++c) count(i + 1, c) = count(i, c);
        const std::size_t k = ev[i].cls;
        const std::size_t from = static_cast<std::size_t>(ev[i].logger) - 1;
        const std::size_t to = (from + 1) % GIBBS_NNODES;
        count(i + 1, from * K + k) -= 1;
        count(i + 1, to * K + k) += 1;
    }

    GibbsStateProbs<T> out;
    out.N.assign(K, 0L);
    for (std::size_t k = 0; k < K; ++k) {
        long m = 0;
        for (std::size_t i = 0; i < total; ++i)
            for (std::size_t nd = 0; nd < GIBBS_NNODES; ++nd)
                if (count(i, nd * K + k) > m) m = count(i, nd * K + k);
        out.N[k] = m;
    }
    for (std::size_t i = 0; i < total; ++i)
        for (std::size_t k = 0; k < K; ++k) count(i, k) += out.N[k];

    // burnin index rationale: see _kb/03-api-layer.md (cpp port notes: infer)
    std::size_t burnin = 0;  // 0-based
    if (data_needed != 0) {
        if (total == data_needed)
            throw InputError(
                "gibbs_analyse_data: the record is exactly data_needed events long, which MATLAB "
                "indexes from zero");
        if (total > data_needed) burnin = total - data_needed - 1;
    }

    // state holding-time aggregation order: see _kb/03-api-layer.md (cpp port notes: infer)
    const T zero = num_traits<T>::from_int(0);
    std::map<std::vector<long>, T> acc;
    for (std::size_t i = burnin; i < total; ++i) {
        std::vector<long> key(K * GIBBS_NNODES);
        for (std::size_t c = 0; c < K * GIBBS_NNODES; ++c) key[c] = count(i, c);
        const T dt = i == 0 ? zero : T(ev[i].t - ev[i - 1].t);
        typename std::map<std::vector<long>, T>::iterator it = acc.find(key);
        if (it == acc.end())
            acc.insert(std::make_pair(key, dt));
        else
            it->second += dt;
    }

    const T obs_length = ev[total - 1].t - ev[burnin].t;
    if (obs_length <= zero) throw NumericError("gibbs_analyse_data: the record has no duration");

    out.states = Matrix<long>(acc.size(), K * GIBBS_NNODES, 0L);
    out.prob.assign(acc.size(), zero);
    std::size_t r = 0;
    for (typename std::map<std::vector<long>, T>::const_iterator it = acc.begin();
         it != acc.end(); ++it, ++r) {
        for (std::size_t c = 0; c < K * GIBBS_NNODES; ++c) out.states(r, c) = it->first[c];
        out.prob[r] = it->second / obs_length;
    }

    out.N0.assign(K, zero);
    for (std::size_t k = 0; k < K; ++k) {
        T s = zero;
        for (std::size_t i = 0; i < out.prob.size(); ++i)
            s += out.prob[i] * num_traits<T>::from_int(out.states(i, K + k));
        out.N0[k] = s;
    }
    return out;
}

/**
 * The deterministic half of one coordinate update: the log normalizing
 * constant along the grid and the normalized slice it implies.
 *
 * @param think_time (K) think times of the delay node
 * @param theta      (K) current demands; theta[index] must lie ON the grid
 * @param testset    (nt x K*GIBBS_NNODES) sampled states
 * @param index      coordinate being updated
 * @param N          (K) populations
 * @param logG_init  running log G at the current theta
 * @param interval   grid step
 * @param range_size upper end of the grid
 */
template <class T>
GibbsSlice<T> gibbs_slice(const std::vector<T>& think_time, const std::vector<T>& theta,
                          const Matrix<long>& testset, std::size_t index,
                          const std::vector<T>& N, const T& logG_init, double interval,
                          const T& range_size) {
    static_assert(num_traits<T>::has_transcendental,
                  "gibbs_slice requires transcendental arithmetic: it integrates "
                  "d log G / d theta along a grid in logarithms");
    using std::exp;
    using std::log;

    const std::size_t K = theta.size();
    if (think_time.size() != K) throw InputError("gibbs_slice: think_time has the wrong length");
    if (N.size() != K) throw InputError("gibbs_slice: N has the wrong length");
    if (index >= K) throw InputError("gibbs_slice: coordinate out of range");
    if (testset.cols() != K * GIBBS_NNODES)
        throw InputError("gibbs_slice: testset has the wrong width");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T eps = detail::gibbs_eps<T>();
    const T step = num_traits<T>::from_double(interval);

    GibbsSlice<T> out;
    const std::size_t n =
        detail::colon_count(num_traits<T>::to_double(range_size), interval);
    if (n == 0) throw NumericError("gibbs_slice: the candidate grid collapsed to nothing");
    out.grid.assign(n, zero);
    for (std::size_t i = 0; i < n; ++i)
        out.grid[i] = num_traits<T>::from_int(static_cast<long>(i)) * step;
    out.logG.assign(n, zero);

    // grid-miss rationale: see _kb/03-api-layer.md (cpp port notes: infer)
    std::size_t ip = 0;
    bool found = false;
    for (std::size_t i = 0; i < n && !found; ++i)
        if (out.grid[i] == theta[index]) {
            ip = i;
            found = true;
        }

    if (found) {
        Matrix<T> L(1, K, zero);
        for (std::size_t j = 0; j < K; ++j) L(0, j) = theta[j];
        out.logG[ip] = logG_init;

        pfqn::AmvaResult<T> res =
            pfqn::pfqn_bs(L, N, think_time, std::vector<pfqn::AmvaSched>());
        Matrix<T> QN = res.QN;
        for (std::size_t i = ip; i-- > 0;) {
            L(0, index) = out.grid[i + 1];
            res = pfqn::pfqn_bs(L, N, think_time, std::vector<pfqn::AmvaSched>(), interval, 1000,
                                QN);
            QN = res.QN;
            const T d = one - QN(0, index) / (out.grid[i + 1] + eps) * step;
            out.logG[i] = d < zero ? out.logG[i + 1] : T(out.logG[i + 1] + log(d));
        }

        L(0, index) = theta[index];
        res = pfqn::pfqn_bs(L, N, think_time, std::vector<pfqn::AmvaSched>());
        QN = res.QN;
        for (std::size_t i = ip + 1; i < n; ++i) {
            L(0, index) = out.grid[i - 1];
            res = pfqn::pfqn_bs(L, N, think_time, std::vector<pfqn::AmvaSched>(), interval, 1000,
                                QN);
            QN = res.QN;
            const T d = one + QN(0, index) / (out.grid[i - 1] + eps) * step;
            out.logG[i] = d < zero ? out.logG[i - 1] : T(out.logG[i - 1] + log(d));
        }
    }

    // Slice log density. The coefficient is the total number of class-index
    // jobs at the queue over the whole test set.
    T coeff = zero;
    for (std::size_t i = 0; i < testset.rows(); ++i)
        coeff += num_traits<T>::from_int(testset(i, K + index));
    if (coeff == zero)
        throw NumericError(
            "gibbs_slice: no sampled state has a job of this class at the queue, so the slice "
            "density is 0 log 0 at the first grid point and carries no information");
    const T nt = num_traits<T>::from_int(static_cast<long>(testset.rows()));

    std::vector<T> lp(n, zero);
    for (std::size_t i = 0; i < n; ++i) lp[i] = coeff * log(out.grid[i]) - out.logG[i] * nt;
    T mx = lp[0];
    for (std::size_t i = 1; i < n; ++i)
        if (lp[i] > mx) mx = lp[i];
    out.prob.assign(n, zero);
    T tot = zero;
    for (std::size_t i = 0; i < n; ++i) {
        out.prob[i] = exp(T(lp[i] - mx));
        tot += out.prob[i];
    }
    if (!(tot > zero)) throw NumericError("gibbs_slice: the slice has no mass");
    for (std::size_t i = 0; i < n; ++i) out.prob[i] = out.prob[i] / tot;

    // empty-find slice-mass rationale: see _kb/03-api-layer.md (cpp port notes: infer)
    T cum = zero;
    const T closed = one - num_traits<T>::from_double(1e-10);
    out.range_size_dim = out.grid[n - 1] * num_traits<T>::from_int(2);
    for (std::size_t i = 0; i < n; ++i) {
        cum += out.prob[i];
        if (cum > closed) {
            out.range_size_dim = out.grid[i] * num_traits<T>::from_int(2);
            break;
        }
    }
    return out;
}

/**
 * Estimated per-class mean demands.
 *
 * @param data    per-class traces
 * @param nbCores number of processors of the queue node
 * @param opts    sample budgets and tolerance
 * @param rng     generator, advanced by the call
 * @return        (K) estimated demands
 */
template <class T>
std::vector<T> infer_gibbs(const std::vector<GibbsTrace<T>>& data, const T& nbCores,
                           const GibbsOptions& opts, pfqn::McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "infer_gibbs requires transcendental arithmetic: it samples a slice of a "
                  "log density built from an iteratively integrated normalizing constant");
    using std::log;

    const std::size_t K = data.size();
    if (K == 0) throw InputError("infer_gibbs: no classes");
    if (opts.block == 0) throw InputError("infer_gibbs: the convergence block must be positive");
    if (opts.nsamples == 0) throw InputError("infer_gibbs: the chain must have samples");
    if (!(opts.tol > 0.0)) throw InputError("infer_gibbs: the tolerance must be positive");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const GibbsStateProbs<T> st = gibbs_analyse_data(data, opts.data_needed);
    const std::size_t ns = st.prob.size();

    // usedCores rationale: see _kb/03-api-layer.md (cpp port notes: infer)
    T used = zero;
    for (std::size_t i = 0; i < ns; ++i) {
        T q = zero;
        for (std::size_t k = 0; k < K; ++k) q += num_traits<T>::from_int(st.states(i, K + k));
        used += (q > nbCores ? nbCores : q) * st.prob[i];
    }
    const T denom = one - st.prob[ns - 1];
    if (denom == zero) throw NumericError("infer_gibbs: the record holds a single state");
    used = used / denom;

    std::vector<T> Nt(K, zero), think_time(K, zero);
    for (std::size_t k = 0; k < K; ++k) {
        Nt[k] = num_traits<T>::from_int(st.N[k]);
        if (data[k].think_obs.empty())
            throw InputError("infer_gibbs: a class has no think-time observations");
        T s = zero;
        for (std::size_t i = 0; i < data[k].think_obs.size(); ++i) s += data[k].think_obs[i];
        const T m = s / num_traits<T>::from_int(static_cast<long>(data[k].think_obs.size()));
        if (m == zero) throw NumericError("infer_gibbs: a class has a zero think-time mean");
        think_time[k] = (Nt[k] - st.N0[k]) / m;
        if (!(think_time[k] > zero))
            throw NumericError("infer_gibbs: a class has a non-positive think time");
    }

    // Test set: states drawn from the empirical distribution by inversion.
    std::vector<T> cum(ns, zero);
    T c = zero;
    for (std::size_t i = 0; i < ns; ++i) {
        c += st.prob[i];
        cum[i] = c;
    }
    Matrix<long> testset(opts.likelihood_sample, K * GIBBS_NNODES, 0L);
    for (std::size_t s = 0; s < opts.likelihood_sample; ++s) {
        const T u = pfqn::mc_uniform<T>(rng);
        std::size_t pick = ns;  // MATLAB indexes an empty find and raises
        for (std::size_t i = 0; i < ns; ++i)
            if (u < cum[i]) {
                pick = i;
                break;
            }
        if (pick == ns)
            throw NumericError("infer_gibbs: the empirical distribution does not sum to one");
        for (std::size_t j = 0; j < K * GIBBS_NNODES; ++j) testset(s, j) = st.states(pick, j);
    }

    // log G of the demand-free model: all jobs at the delay node.
    T logG = zero;
    for (std::size_t k = 0; k < K; ++k) {
        logG += Nt[k] * log(think_time[k]);
        for (long j = 1; j <= st.N[k]; ++j) logG -= log(num_traits<T>::from_int(j));
    }

    std::vector<T> range_size(K, one);
    std::vector<T> theta(K, zero);
    Matrix<T> smpl(opts.nsamples, K, zero);
    std::vector<T> demand_old(K, zero);
    const std::size_t nblocks = static_cast<std::size_t>(
        std::floor(static_cast<double>(opts.nsamples) / static_cast<double>(opts.block) + 0.5));
    std::size_t si = 0;  // number of samples drawn so far

    for (std::size_t b = 1; b <= nblocks; ++b) {
        for (std::size_t s = 0; s < opts.block && si < opts.nsamples; ++s) {
            for (std::size_t h = 0; h < K; ++h) {
                // theta: this sweep's coordinates below h, the previous
                // sweep's at and above h (zeros on the first sweep).
                for (std::size_t j = 0; j < h; ++j) theta[j] = smpl(si, j);
                for (std::size_t j = h; j < K; ++j)
                    theta[j] = si == 0 ? zero : smpl(si - 1, j);

                const GibbsSlice<T> sl =
                    gibbs_slice(think_time, theta, testset, h, Nt, logG, opts.tol, range_size[h]);

                const T u = pfqn::mc_uniform<T>(rng);
                T cc = zero;
                std::size_t pick = sl.grid.size();
                for (std::size_t i = 0; i < sl.grid.size(); ++i) {
                    cc += sl.prob[i];
                    if (u < cc) {
                        pick = i;
                        break;
                    }
                }
                if (pick == sl.grid.size()) {
                    // MATLAB's empty-find branch: keep the current value and
                    // the running constant untouched.
                    smpl(si, h) = theta[h];
                } else {
                    smpl(si, h) = sl.grid[pick];
                    logG = sl.logG[pick];
                }
                // grid-width doubling rationale: see _kb/03-api-layer.md (cpp port notes: infer)
                range_size[h] = sl.range_size_dim * num_traits<T>::from_int(2);
            }
            ++si;
        }

        if (b == 2) {
            for (std::size_t k = 0; k < K; ++k) {
                T s = zero;
                for (std::size_t i = opts.block; i < si; ++i) s += smpl(i, k);
                demand_old[k] = s / num_traits<T>::from_int(static_cast<long>(si - opts.block));
            }
        } else if (b > 2) {
            const std::size_t lo = (b - 1) * opts.block;
            std::vector<T> demand_now(K, zero);
            const T nb = num_traits<T>::from_int(static_cast<long>(si - lo));
            const T bp1 = num_traits<T>::from_int(static_cast<long>(b + 1));
            const T bb = num_traits<T>::from_int(static_cast<long>(b));
            for (std::size_t k = 0; k < K; ++k) {
                T s = zero;
                for (std::size_t i = lo; i < si; ++i) s += smpl(i, k);
                const T m = s / nb;
                demand_now[k] = m / bp1 + demand_old[k] / bp1 * bb;
            }
            T rel = zero;
            bool ok = true;
            for (std::size_t k = 0; k < K; ++k) {
                if (demand_old[k] == zero) {
                    ok = false;
                    break;
                }
                rel += num_abs(T((demand_now[k] - demand_old[k]) / demand_old[k]));
            }
            if (ok) rel = rel / num_traits<T>::from_int(static_cast<long>(K));
            if (ok && num_traits<T>::to_double(rel) < opts.tol) break;
            demand_old = demand_now;
        }
    }

    // MATLAB averages the second half of the chain, dropping the last sample.
    const std::size_t nb = si == 0 ? 0 : si - 1;
    const std::size_t lo = static_cast<std::size_t>(
        std::floor(static_cast<double>(nb) / 2.0 + 0.5));  // round(nb/2)+1, 0-based
    if (nb == 0 || lo >= nb) throw NumericError("infer_gibbs: the chain produced no usable tail");
    std::vector<T> demand(K, zero);
    for (std::size_t k = 0; k < K; ++k) {
        T s = zero;
        for (std::size_t i = lo; i < nb; ++i) s += smpl(i, k) * used;
        demand[k] = s / num_traits<T>::from_int(static_cast<long>(nb - lo));
    }
    return demand;
}

/** MATLAB's three-argument form, with its hard-coded budgets. */
template <class T>
std::vector<T> infer_gibbs(const std::vector<GibbsTrace<T>>& data, const T& nbCores, double tol,
                           pfqn::McRng& rng) {
    GibbsOptions o;
    o.tol = tol;
    return infer_gibbs(data, nbCores, o, rng);
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_GIBBS_H
