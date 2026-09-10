/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_QLEN_JOINT_MOMENTS_H
#define LINE_API_PFQN_QLEN_JOINT_MOMENTS_H

/**
 * Joint moments of the queue-length vector of a closed product-form network,
 * obtained from normalizing constants.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_qlen_joint_moments.m,
 * cross-checked against jar/src/main/java/jline/api/pfqn/nc/
 * Pfqn_qlen_joint_moments.java.
 *
 * The coordinates are (station, class) pairs. Two pairs sharing a class give the
 * cross-station covariance of that class; two pairs sharing a station give the
 * cross-class covariance at that station, which is what a class-oriented method
 * of moments (pfqn_comomrm and its relatives) is positioned to deliver. Two
 * exact routes reach the joint survival array, and both end in the same
 * conversion, the tail edge of the house of moments (api/moment) followed by the
 * joint central-moment and cumulant conversions:
 *
 *   - SINGLE CLASS (R = 1), route 'tail'. The survival probabilities are ratios
 *     of normalizing constants of the network itself,
 *
 *       P(n_i >= k_i for all i) = (prod_i L_i^k_i) G(N - sum_i k_i) / G(N)
 *
 *     which holds because a load-independent single-class station has the
 *     geometric occupancy L_i^n. Only N+1 constants of the ORIGINAL model are
 *     needed, which is why any normalizing-constant algorithm serves it.
 *
 *   - MULTICLASS, route 'pmf'. The geometric factorization fails, since a
 *     multiclass load-independent station carries the multinomial occupancy
 *     f_i(n_i) = |n_i|! prod_r L_ir^n_ir / n_ir!. What holds instead is the joint
 *     law of the selected stations in terms of the COMPLEMENTARY network, the
 *     model with those stations deleted and the think times kept,
 *
 *       P(n_i = m_i, i in S) = prod_i f_i(m_i) G_(S^c)(N - sum_i m_i) / G(N).
 *
 *     The survival array is the reverse cumulative sum of that array, exactly,
 *     since the box covers the support.
 *
 * Neither the factorial nor the raw moments have a one-constant closed form; the
 * survival array is the queue-length functional that does. The
 * normalizing-constant algorithm is INJECTED rather than called at a fixed site:
 * the whole set of populations is known before any evaluation, so it is emitted
 * in one batch and an algorithm that produces several constants in one pass
 * serves it without recomputation.
 *
 * THE DEFAULT METHOD IS EXACT, NOT ADAPTIVE. The point of this routine is an
 * exact moment array; an approximate normalizing constant would silently make
 * every moment approximate. The reference builds its options from
 * SolverNC.defaultOptions and overrides only the method, which this port
 * reproduces by defaulting NcMethod to Exact.
 *
 * Reference: M. Reiser and S. S. Lavenberg, "Mean-value analysis of closed
 * multichain queuing networks", JACM 27(2):313-322, 1980.
 *
 * Arithmetic: TRANSCENDENTAL. The survival array is assembled in the log domain
 * so that no ratio of normalizing constants overflows.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/api/moment/moment_joint.h"
#include "line/api/moment/moment_tensor.h"
#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_nc.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Which survival-array identity is used. */
enum class QlenJointRoute {
    Auto,  ///< tail for a single class, pmf otherwise
    Tail,  ///< the geometric single-class identity
    Pmf    ///< the complementary-network joint law
};

/** Everything pfqn_qlen_joint_moments reports. */
template <class T>
struct QlenJointMomentsResult {
    moment::MomentTensor<T> tail;
    moment::MomentTensor<T> binomial;
    moment::MomentTensor<T> factorial;
    moment::MomentTensor<T> raw;
    moment::MomentTensor<T> central;
    moment::MomentTensor<T> cumulant;
    std::vector<T> mean;  ///< (P) mean of each coordinate
    Matrix<T> cov;        ///< (P x P) covariance of the coordinates
    // info
    std::string route;
    std::size_t points = 0;  ///< distinct populations requested
    std::size_t served = 0;  ///< how many the injected source supplied
    std::size_t evals = 0;   ///< how many needed a pfqn_nc call
    bool exact = true;
    std::vector<std::pair<std::size_t, std::size_t>> pairs;  ///< 0-based (station, class)
    std::vector<std::size_t> dims;
};

/**
 * Injected source of log G. Invoked ONCE per network as lGsrc(Lsub, pops), pops
 * being a list of populations; it must return one value per population and may
 * return NaN where it cannot serve, which is then filled in by pfqn_nc. The
 * 'pmf' route queries the COMPLEMENTARY network, so a source must answer for
 * whichever demand matrix it is handed.
 */
template <class T>
using QlenJointLgSource =
    std::function<std::vector<double>(const Matrix<T>&, const std::vector<std::vector<int>>&)>;

namespace detail {

/** Advance a 0-based multi-index, first dimension fastest. */
inline bool qlen_odometer(std::vector<std::size_t>& a, const std::vector<std::size_t>& dims) {
    for (std::size_t l = 0; l < dims.size(); ++l) {
        if (++a[l] < dims[l]) return true;
        a[l] = 0;
    }
    return false;
}

/**
 * MomentTensor drops trailing singleton dimensions, so a d-length multi-index
 * has to be truncated to the tensor order; the dropped entries are always zero
 * because their extent is one.
 */
template <class T>
std::vector<std::size_t> qlen_trunc(const moment::MomentTensor<T>& A,
                                    const std::vector<std::size_t>& ord) {
    std::vector<std::size_t> out(ord.begin(), ord.begin() + std::min(ord.size(), A.order()));
    out.resize(A.order(), 0);
    return out;
}

/** Index of a population vector in the deduplicated request list. */
inline std::size_t qlen_findrow(const std::vector<std::vector<int>>& rows,
                                const std::vector<int>& key) {
    for (std::size_t p = 0; p < rows.size(); ++p)
        if (rows[p] == key) return p;
    throw NumericError("pfqn_qlen_joint_moments: a requested population was not evaluated");
}

/** Log normalizing constant of a pure-delay network, prod_r Z_r^n_r / n_r!. */
template <class T>
double qlen_delay_lg(const std::vector<T>& Z, const std::vector<int>& n) {
    using std::log;
    double lg = 0.0;
    for (std::size_t r = 0; r < n.size(); ++r) {
        if (n[r] == 0) continue;
        if (Z[r] <= num_traits<T>::from_int(0)) return -std::numeric_limits<double>::infinity();
        lg += static_cast<double>(n[r]) * num_traits<T>::log_as_double(Z[r]) -
              num_traits<T>::to_double(
                  num_lgamma<T>(num_traits<T>::from_int(static_cast<long>(n[r]) + 1)));
    }
    return lg;
}

/** Evaluate log G at a batch of populations, honouring the injected source. */
template <class T>
std::vector<double> qlen_batch_lg(const Matrix<T>& Lsub, const std::vector<std::vector<int>>& pops,
                                  const Matrix<T>& Z, const QlenJointLgSource<T>& lGsrc,
                                  NcMethod method, const NcOptions& nopt, std::size_t& served,
                                  std::size_t& evals) {
    const std::size_t P = pops.size();
    std::vector<double> lg(P, std::numeric_limits<double>::quiet_NaN());
    served = 0;
    if (lGsrc) {
        lg = lGsrc(Lsub, pops);
        if (lg.size() != P)
            throw InputError(
                "pfqn_qlen_joint_moments: the lGsrc source must return one value per requested "
                "population");
        for (std::size_t p = 0; p < P; ++p)
            if (std::isfinite(lg[p])) ++served;
    }
    evals = 0;
    for (std::size_t p = 0; p < P; ++p) {
        if (std::isfinite(lg[p])) continue;
        const std::vector<T> lambda(pops[p].size(), num_traits<T>::from_int(0));
        lg[p] = pfqn_nc(lambda, Lsub, pops[p], Z, method, num_traits<T>::from_int(0), nopt).lG;
        ++evals;
    }
    return lg;
}

}  // namespace detail

/**
 * @param L      (M x R) demands of the QUEUEING stations; delay stations belong
 *               in Z, their marginals following a different law
 * @param N      (R) population vector
 * @param Z      (R) think times, zeros for none
 * @param pairs  0-based (station, class) coordinates, one per dimension of the
 *               returned arrays; empty for every class of every station
 * @param route  which survival identity to use
 * @param lGsrc  injected source of log G; empty to call pfqn_nc throughout
 * @param method the normalizing-constant algorithm handed to pfqn_nc
 * @param nopt   sample count and seed the estimators read
 */
template <class T>
QlenJointMomentsResult<T> pfqn_qlen_joint_moments(
    const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
    const std::vector<std::pair<std::size_t, std::size_t>>& pairs, QlenJointRoute route,
    const QlenJointLgSource<T>& lGsrc, NcMethod method, const NcOptions& nopt) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_qlen_joint_moments assembles the survival array in the log domain and "
                  "needs transcendental arithmetic");
    using std::exp;
    using std::log;
    const std::size_t M = L.rows();
    const std::size_t R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    if (N.size() != R)
        throw InputError(
            "pfqn_qlen_joint_moments: the population vector N must have one entry per class");
    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(R, zero);
    if (Zv.size() != R)
        throw InputError(
            "pfqn_qlen_joint_moments: the think time vector Z must have one entry per class");
    Matrix<T> Zm(1, R, zero);
    for (std::size_t r = 0; r < R; ++r) Zm(0, r) = Zv[r];

    std::vector<std::pair<std::size_t, std::size_t>> pr = pairs;
    if (pr.empty()) {
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) pr.push_back(std::make_pair(i, r));
    }
    for (std::size_t j = 0; j < pr.size(); ++j)
        if (pr[j].first >= M || pr[j].second >= R)
            throw InputError("pfqn_qlen_joint_moments: a (station,class) pair is out of range");
    {
        std::vector<std::pair<std::size_t, std::size_t>> s = pr;
        std::sort(s.begin(), s.end());
        if (std::unique(s.begin(), s.end()) != s.end())
            throw InputError("pfqn_qlen_joint_moments: the (station,class) pairs must be distinct");
    }

    if (route == QlenJointRoute::Auto) route = (R == 1) ? QlenJointRoute::Tail : QlenJointRoute::Pmf;
    if (route == QlenJointRoute::Tail && R > 1)
        throw InputError(
            "pfqn_qlen_joint_moments: the tail route needs the geometric occupancy of a "
            "single-class load-independent station; with several classes the multinomial factor "
            "breaks the survival identity, so use the pmf route");

    const std::size_t d = pr.size();
    std::vector<std::size_t> dims(d);
    for (std::size_t j = 0; j < d; ++j)
        dims[j] = static_cast<std::size_t>(N[pr[j].second]) + 1;

    QlenJointMomentsResult<T> res;
    res.pairs = pr;
    res.dims = dims;
    res.route = (route == QlenJointRoute::Tail) ? "tail" : "pmf";

    moment::MomentTensor<T> tail(dims);
    std::size_t served = 0, evals = 0, points = 0;

    if (route == QlenJointRoute::Tail) {
        // the whole population set is known up front: N minus the total order
        std::vector<std::vector<int>> need;
        std::vector<std::size_t> a(d, 0);
        do {
            long s = 0;
            for (std::size_t j = 0; j < d; ++j) s += static_cast<long>(a[j]);
            if (static_cast<long>(N[0]) - s >= 0) {
                std::vector<int> row(1, static_cast<int>(static_cast<long>(N[0]) - s));
                if (std::find(need.begin(), need.end(), row) == need.end()) need.push_back(row);
            }
        } while (detail::qlen_odometer(a, dims));
        if (std::find(need.begin(), need.end(), N) == need.end()) need.push_back(N);
        std::sort(need.begin(), need.end());

        std::size_t sv = 0, ev = 0;
        const std::vector<double> lg =
            detail::qlen_batch_lg(L, need, Zm, lGsrc, method, nopt, sv, ev);
        served = sv;
        evals = ev;
        points = need.size();
        const double lgN = lg[detail::qlen_findrow(need, N)];

        std::fill(a.begin(), a.end(), static_cast<std::size_t>(0));
        std::size_t ia = 0;
        do {
            long s = 0;
            for (std::size_t j = 0; j < d; ++j) s += static_cast<long>(a[j]);
            if (static_cast<long>(N[0]) - s >= 0) {
                double acc = 0.0;
                bool ok = true;
                for (std::size_t j = 0; j < d; ++j) {
                    if (a[j] == 0) continue;
                    const T Lij = L(pr[j].first, pr[j].second);
                    if (Lij <= zero) {
                        ok = false;
                        break;
                    }
                    acc += static_cast<double>(a[j]) * num_traits<T>::log_as_double(Lij);
                }
                if (ok) {
                    std::vector<int> key(1, static_cast<int>(static_cast<long>(N[0]) - s));
                    tail.data[ia] = num_traits<T>::from_double(
                        std::exp(acc + lg[detail::qlen_findrow(need, key)] - lgN));
                }
            }
            ++ia;
        } while (detail::qlen_odometer(a, dims));
    } else {
        // the joint law of the selected stations needs every class of those
        // stations, so the internal box runs over (station,class) and the
        // requested pairs are marginalized out of it afterwards
        std::vector<std::size_t> stations;
        for (std::size_t j = 0; j < d; ++j) stations.push_back(pr[j].first);
        std::sort(stations.begin(), stations.end());
        stations.erase(std::unique(stations.begin(), stations.end()), stations.end());
        std::vector<std::pair<std::size_t, std::size_t>> coords;
        for (std::size_t si = 0; si < stations.size(); ++si)
            for (std::size_t r = 0; r < R; ++r) coords.push_back(std::make_pair(stations[si], r));
        const std::size_t dc = coords.size();
        std::vector<std::size_t> cdims(dc);
        for (std::size_t j = 0; j < dc; ++j)
            cdims[j] = static_cast<std::size_t>(N[coords[j].second]) + 1;

        std::vector<std::size_t> keep;
        for (std::size_t i = 0; i < M; ++i)
            if (std::find(stations.begin(), stations.end(), i) == stations.end()) keep.push_back(i);
        Matrix<T> Lsub(keep.size(), R, zero);
        for (std::size_t k = 0; k < keep.size(); ++k)
            for (std::size_t r = 0; r < R; ++r) Lsub(k, r) = L(keep[k], r);

        std::vector<std::vector<int>> need;
        std::vector<std::size_t> a(dc, 0);
        do {
            std::vector<int> n = N;
            for (std::size_t j = 0; j < dc; ++j) n[coords[j].second] -= static_cast<int>(a[j]);
            bool ok = true;
            for (std::size_t r = 0; r < R; ++r)
                if (n[r] < 0) ok = false;
            if (ok && std::find(need.begin(), need.end(), n) == need.end()) need.push_back(n);
        } while (detail::qlen_odometer(a, cdims));
        std::sort(need.begin(), need.end());

        std::vector<double> lgc(need.size(), 0.0);
        if (keep.empty()) {
            for (std::size_t p = 0; p < need.size(); ++p) lgc[p] = detail::qlen_delay_lg(Zv, need[p]);
            served = need.size();
            evals = 0;
        } else {
            std::size_t sv = 0, ev = 0;
            lgc = detail::qlen_batch_lg(Lsub, need, Zm, lGsrc, method, nopt, sv, ev);
            served = sv;
            evals = ev;
        }
        points = need.size();

        std::size_t sv0 = 0, ev0 = 0;
        const QlenJointLgSource<T> none;
        const std::vector<std::vector<int>> onlyN(1, N);
        const std::vector<double> lgNv =
            detail::qlen_batch_lg(L, onlyN, Zm, none, method, nopt, sv0, ev0);
        const double lgN = lgNv[0];
        evals += ev0;

        moment::MomentTensor<T> marg(dims);
        std::fill(a.begin(), a.end(), static_cast<std::size_t>(0));
        do {
            std::vector<int> n = N;
            for (std::size_t j = 0; j < dc; ++j) n[coords[j].second] -= static_cast<int>(a[j]);
            bool inbox = true;
            for (std::size_t r = 0; r < R; ++r)
                if (n[r] < 0) inbox = false;
            if (!inbox) continue;
            const double gc = lgc[detail::qlen_findrow(need, n)];
            if (!std::isfinite(gc)) continue;
            double acc = gc - lgN;
            bool ok = true;
            for (std::size_t si = 0; si < stations.size() && ok; ++si) {
                const std::size_t i = stations[si];
                long tot = 0;
                for (std::size_t j = 0; j < dc; ++j)
                    if (coords[j].first == i) tot += static_cast<long>(a[j]);
                acc += num_traits<T>::to_double(
                    detail::num_lgamma<T>(num_traits<T>::from_int(tot + 1)));
                for (std::size_t j = 0; j < dc; ++j) {
                    if (coords[j].first != i || a[j] == 0) continue;
                    const T Lir = L(i, coords[j].second);
                    if (Lir <= zero) {
                        ok = false;
                        break;
                    }
                    acc += static_cast<double>(a[j]) * num_traits<T>::log_as_double(Lir) -
                           num_traits<T>::to_double(detail::num_lgamma<T>(
                               num_traits<T>::from_int(static_cast<long>(a[j]) + 1)));
                }
            }
            if (!ok) continue;
            std::vector<std::size_t> sub(d, 0);
            for (std::size_t j = 0; j < d; ++j)
                for (std::size_t jc = 0; jc < dc; ++jc)
                    if (coords[jc] == pr[j]) sub[j] = a[jc];
            marg.at(detail::qlen_trunc(marg, sub)) += num_traits<T>::from_double(std::exp(acc));
        } while (detail::qlen_odometer(a, cdims));

        // the survival array is the reverse cumulative sum along every mode
        tail = marg;
        for (std::size_t mode = 0; mode < tail.order(); ++mode) {
            const std::size_t n = tail.sz[mode];
            std::size_t stride = 1;
            for (std::size_t l = 0; l < mode; ++l) stride *= tail.sz[l];
            const std::size_t outer = tail.numel() / (n * stride);
            for (std::size_t o = 0; o < outer; ++o)
                for (std::size_t s = 0; s < stride; ++s) {
                    const std::size_t base = o * n * stride + s;
                    for (std::size_t k = n - 1; k-- > 0;)
                        tail.data[base + k * stride] += tail.data[base + (k + 1) * stride];
                }
        }
    }

    res.tail = tail;
    res.binomial = moment::moment_joint_binomial_from_tail(tail);
    res.factorial = moment::moment_joint_factorial_from_binomial(res.binomial);
    res.raw = moment::moment_joint_raw_from_factorial(res.factorial);
    res.central = moment::moment_joint_central_from_raw(res.raw);
    res.cumulant = moment::moment_joint_cumulant_from_raw(res.raw);

    res.mean.assign(d, zero);
    res.cov = Matrix<T>(d, d, zero);
    for (std::size_t j = 0; j < d; ++j) {
        std::vector<std::size_t> e(d, 0);
        e[j] = 1;
        res.mean[j] = res.raw.at(detail::qlen_trunc(res.raw, e));
        for (std::size_t l = 0; l < d; ++l) {
            std::vector<std::size_t> aa(d, 0);
            ++aa[j];
            ++aa[l];
            res.cov(j, l) = res.cumulant.at(detail::qlen_trunc(res.cumulant, aa));
        }
    }
    res.points = points;
    res.served = served;
    res.evals = evals;
    return res;
}

/** MATLAB defaults: every coordinate, the automatic route, no injected source. */
template <class T>
QlenJointMomentsResult<T> pfqn_qlen_joint_moments(const Matrix<T>& L, const std::vector<int>& N,
                                                  const std::vector<T>& Z) {
    return pfqn_qlen_joint_moments(L, N, Z, std::vector<std::pair<std::size_t, std::size_t>>(),
                                   QlenJointRoute::Auto, QlenJointLgSource<T>(), NcMethod::Exact,
                                   NcOptions());
}

/** MATLAB defaults with an explicit coordinate list. */
template <class T>
QlenJointMomentsResult<T> pfqn_qlen_joint_moments(
    const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
    const std::vector<std::pair<std::size_t, std::size_t>>& pairs) {
    return pfqn_qlen_joint_moments(L, N, Z, pairs, QlenJointRoute::Auto, QlenJointLgSource<T>(),
                                   NcMethod::Exact, NcOptions());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_QLEN_JOINT_MOMENTS_H
