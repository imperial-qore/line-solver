/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MVAOI_H
#define LINE_API_PFQN_MVAOI_H

/**
 * Mean-value analysis of a closed network with order-independent (OI) stations,
 * the composition-dependent generalization of Conditional MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mvaoi.m. The network is an
 * aggregated delay node, any number of load-independent single-server queues,
 * and any number of OI / pass-and-swap stations with empty swap graph. It
 * returns the same exact throughputs and queue lengths as pfqn_ncoi WITHOUT
 * forming any normalizing constant or joint marginal.
 *
 * The recursion carries, per OI station i, the shift vector s_i, the OI
 * occupancy already committed at the bottom of that station; S is the K x R
 * matrix of those rows and Nn = N - sum_i s_i the jobs still free. For a single
 * OI station and no LI queue,
 *
 *   Q^{(S)}(Nn)   = sum_r U_r^{(S)}(Nn) ( e_r + Q^{(S + e_r@@i)}(Nn - e_r) )
 *   U_r^{(S)}(Nn) = D_r^{(S)}(Nn) X_r^{(S)}(Nn)
 *   D_r^{(S)}(Nn) = (1/mu_i(s_i + e_r)) rho^{(S)}_{i,r}(Nn - e_r),          Nn_r = 1
 *   D_r^{(S)}(Nn) = [X_r^{(S)}/X_r^{(S+e_r@@i)}](Nn - e_r) D_r^{(S)}(Nn-e_r), Nn_r >= 2
 *   rho^{(S)}_{i,r}(M) = rho^{(S)}_{i,r}(M - e_s) X_s^{(S)}(M)/X_s^{(S+e_r@@i)}(M)
 *
 * with each OI station keeping its own D, rho and Q driven by the common
 * throughput, and the population-conservation identity aggregating every
 * station. States (S, Nn) are processed by increasing sum(Nn), so every
 * reference lands at a strictly smaller free population.
 *
 * THE THROUGHPUT CLOSURE IS NOT A LINEAR SOLVE, and that is deliberate in the
 * reference. The OI queue recurrence couples X_s (s != r) through the
 * off-diagonal of A, so A X = Nn is not a per-class Little's-law ratio as in
 * canonical CMVA. Rather than solving the system, the reference decouples it
 * with the product-form throughput-ratio identity at fixed shift,
 * X_s(Nn)/X_r(Nn) = X_s(Nn-e_r)/X_r(Nn-e_s), giving a scalar per-class formula
 * fed by one-job-less throughputs already in the cache. The port keeps that
 * form; substituting a linear solve would change the numbers.
 *
 * Soi, the mean number of IN-SERVICE jobs per class, is the only output that is
 * not a pure mean-value quantity. It is a distributional statistic and comes
 * from the OI count marginal
 *   pM_i(n|k) = (1/mu_i(n)) sum_r X_r(k) pM_i(n - e_r | k - e_r),
 *   pM_i(0|k) = 1 - sum_{n != 0} pM_i(n|k),
 * assembled from the zero-shift throughputs the mean-value recursion has
 * already cached, weighted by pfqn_oi_insvc's E[sir_r | n]. Still no
 * normalizing constant. It is computed only when `want_soi` is set, mirroring
 * the reference's nargout >= 5 guard.
 *
 * Arithmetic: EXACT-CAPABLE. Additions, multiplications, divisions and
 * comparisons in the field of the inputs; no logarithm, no tolerance, no
 * iteration to convergence. The OI rate handle must return a T.
 *
 * REFERENCE DEFECTS: none found. Agreement with pfqn_mvaoi_marg, which reaches
 * the same numbers by an entirely different (marginal-distribution, iterated)
 * route, is to ~2e-14 on every model tried, the gap being that routine's
 * fixed-point tolerance rather than a disagreement.
 */

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <vector>

#include "line/api/pfqn/pfqn_comb_common.h"
#include "line/api/pfqn/pfqn_oi_insvc.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_mvaoi, mirroring [X, Qoi, Qli, Qdelay, Soi]. */
template <class T>
struct MvaoiResult {
    std::vector<T> X;       ///< (R) per-class throughput
    Matrix<T> Qoi;          ///< (K x R) per-class queue length at each OI station
    Matrix<T> Qli;          ///< (J x R) per-class queue length at each LI queue
    std::vector<T> Qdelay;  ///< (R) per-class queue length at the delay node
    Matrix<T> Soi;          ///< (K x R) per-class mean in-service jobs, if requested
};

namespace detail {

/** Lattice key for a state (S, Nn); S is K x R row-major, Nn is R. */
inline std::vector<int> mvaoi_key(const Matrix<int>& S, const std::vector<int>& Nn) {
    std::vector<int> key;
    key.reserve(S.rows() * S.cols() + Nn.size());
    for (std::size_t i = 0; i < S.rows(); ++i)
        for (std::size_t r = 0; r < S.cols(); ++r) key.push_back(S(i, r));
    key.insert(key.end(), Nn.begin(), Nn.end());
    return key;
}

}  // namespace detail

/**
 * @param Z        (R) think-time demands of the aggregated delay node
 * @param N        (R) closed populations
 * @param mu       (K) OI rate handles; mu[i](n) is the total rate of station i
 *                 at the per-class occupancy n
 * @param Dli      (J x R) demands of the load-independent single-server queues
 * @param visits   (K x R) per-OI-station class visit ratios v_{i,r}; they enter the
 *                 class-r demand base case theta_{i,r}(N_r=1) = v_{i,r}/mu_i(...),
 *                 the N_r >= 2 ratio case cancelling them. Empty for unit visits;
 *                 ms-promoted stations pass ones, their visits already folded into
 *                 the rate handle by the caller
 * @param want_soi compute Soi (the reference's nargout >= 5 branch)
 */
template <class T>
MvaoiResult<T> pfqn_mvaoi(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli, const Matrix<T>& visits, bool want_soi) {
    const std::size_t R = N.size();
    if (R == 0) throw InputError("pfqn_mvaoi: empty population vector");
    if (Z.size() != R) throw InputError("pfqn_mvaoi: Z and N disagree on the class count");
    if (mu.empty()) throw InputError("pfqn_mvaoi: mu must be a nonempty list of OI rate handles");
    const std::size_t K = mu.size();
    for (std::size_t i = 0; i < K; ++i)
        if (!mu[i]) throw InputError("pfqn_mvaoi: each mu[i] must be callable");
    const std::size_t J = Dli.empty() ? 0 : Dli.rows();
    if (J > 0 && Dli.cols() != R)
        throw InputError("pfqn_mvaoi: Dli and N disagree on the class count");
    if (!visits.empty() && (visits.rows() != K || visits.cols() != R))
        throw InputError("pfqn_mvaoi: visits must be K x R");
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] < 0) throw InputError("pfqn_mvaoi: pfqn_mvaoi requires finite closed populations");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    std::map<std::vector<int>, std::vector<T>> Xc;   // X^{(S)}(Nn)
    std::map<std::vector<int>, Matrix<T>> Qlc;       // Qli^{(S)}(Nn)
    std::vector<std::map<std::vector<int>, std::vector<T>>> Dc(K), Qc(K);
    std::vector<std::map<std::vector<int>, std::vector<std::pair<bool, T>>>> Rc(K);

    const Matrix<int> zeroS(K, R, 0);

    // (S,Nn) composition rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<std::vector<std::vector<int>>> comps(R);
    for (std::size_t r = 0; r < R; ++r)
        comps[r] = multichoose_rows(static_cast<int>(K) + 2, N[r]);
    std::vector<std::pair<int, std::pair<Matrix<int>, std::vector<int>>>> states;
    {
        std::vector<std::size_t> idx(R, 0);
        bool more = true;
        while (more) {
            Matrix<int> S(K, R, 0);
            std::vector<int> Nn(R, 0);
            int tot = 0;
            for (std::size_t r = 0; r < R; ++r) {
                const std::vector<int>& c = comps[r][idx[r]];
                for (std::size_t i = 0; i < K; ++i) S(i, r) = c[i];
                Nn[r] = c[K];
                tot += Nn[r];
            }
            states.push_back(std::make_pair(tot, std::make_pair(S, Nn)));
            std::size_t d = 0;
            for (; d < R; ++d) {
                if (++idx[d] < comps[d].size()) break;
                idx[d] = 0;
            }
            more = (d < R);
        }
    }
    std::stable_sort(states.begin(), states.end(),
                     [](const std::pair<int, std::pair<Matrix<int>, std::vector<int>>>& x,
                        const std::pair<int, std::pair<Matrix<int>, std::vector<int>>>& y) {
                         return x.first < y.first;
                     });

    // rho^{(S)}_{i,r}(M), memoized per OI station
    std::function<T(std::size_t, std::size_t, const Matrix<int>&, const std::vector<int>&)> rho_fn =
        [&](std::size_t i, std::size_t r, const Matrix<int>& S, const std::vector<int>& Mv) -> T {
        const std::vector<int> rkey = detail::mvaoi_key(S, Mv);
        auto it = Rc[i].find(rkey);
        if (it != Rc[i].end() && it->second[r].first) return it->second[r].second;
        int tot = 0;
        for (int x : Mv) tot += x;
        if (tot == 0) {
            if (it == Rc[i].end())
                it = Rc[i].insert(std::make_pair(rkey, std::vector<std::pair<bool, T>>(
                                                          R, std::make_pair(false, zero))))
                         .first;
            it->second[r] = std::make_pair(true, one);
            return one;
        }
        std::size_t rs = R;
        for (std::size_t t = 0; t < R; ++t)
            if (t != r && Mv[t] > 0) {
                rs = t;
                break;
            }
        if (rs == R) throw NumericError("pfqn_mvaoi: rho recursion has no class to decrement");
        const std::vector<T>& xu = Xc.at(detail::mvaoi_key(S, Mv));
        Matrix<int> Sp = S;
        Sp(i, r) += 1;
        const std::vector<T>& xu2 = Xc.at(detail::mvaoi_key(Sp, Mv));
        T ratio = zero;
        if (xu2[rs] > zero) ratio = xu[rs] / xu2[rs];
        std::vector<int> Mm = Mv;
        Mm[rs] -= 1;
        const T v = rho_fn(i, r, S, Mm) * ratio;
        it = Rc[i].find(rkey);
        if (it == Rc[i].end())
            it = Rc[i].insert(std::make_pair(
                                  rkey, std::vector<std::pair<bool, T>>(R, std::make_pair(false, zero))))
                     .first;
        it->second[r] = std::make_pair(true, v);
        return v;
    };

    for (std::size_t p = 0; p < states.size(); ++p) {
        const Matrix<int>& S = states[p].second.first;
        const std::vector<int>& Nn = states[p].second.second;
        const std::vector<int> key = detail::mvaoi_key(S, Nn);
        if (states[p].first == 0) {
            // No free jobs. Stored for EVERY shift S, not only S = 0, because
            // Qsub references (S + e_s@@i, Nn - e_s) and reaches these.
            Xc[key] = std::vector<T>(R, zero);
            Qlc[key] = Matrix<T>(J, R, zero);
            for (std::size_t i = 0; i < K; ++i) {
                Dc[i][key] = std::vector<T>(R, zero);
                Qc[i][key] = std::vector<T>(R, zero);
            }
            continue;
        }

        Matrix<T> Dt(K, R, zero);
        std::vector<Matrix<T>> Qsub(K, Matrix<T>(R, R, zero));
        for (std::size_t i = 0; i < K; ++i) {
            for (std::size_t r = 0; r < R; ++r) {
                if (Nn[r] == 0) continue;
                std::vector<int> Nr = Nn;
                Nr[r] -= 1;
                if (Nn[r] == 1) {
                    std::vector<int> occ(R, 0);
                    for (std::size_t t = 0; t < R; ++t) occ[t] = S(i, t);
                    occ[r] += 1;
                    const T mur = mu[i](occ);
                    const T vir = visits.empty() ? one : visits(i, r);
                    if (mur > zero) Dt(i, r) = (vir / mur) * rho_fn(i, r, S, Nr);
                } else {
                    Matrix<int> Sp = S;
                    Sp(i, r) += 1;
                    const std::vector<T>& xs = Xc.at(detail::mvaoi_key(S, Nr));
                    const std::vector<T>& xs2 = Xc.at(detail::mvaoi_key(Sp, Nr));
                    if (xs2[r] > zero) {
                        const std::vector<T>& Dprev = Dc[i].at(detail::mvaoi_key(S, Nr));
                        Dt(i, r) = (xs[r] / xs2[r]) * Dprev[r];
                    }
                }
            }
            for (std::size_t s = 0; s < R; ++s) {
                if (Nn[s] == 0) continue;
                Matrix<int> Ss = S;
                Ss(i, s) += 1;
                std::vector<int> Ns = Nn;
                Ns[s] -= 1;
                const std::vector<T>& q = Qc[i].at(detail::mvaoi_key(Ss, Ns));
                for (std::size_t t = 0; t < R; ++t) Qsub[i](s, t) = q[t];
            }
        }

        // LI-queue arrival-theorem coefficients
        Matrix<T> betaLI(J, R, zero);
        for (std::size_t r = 0; r < R; ++r) {
            if (Nn[r] == 0) continue;
            std::vector<int> Nr = Nn;
            Nr[r] -= 1;
            const Matrix<T>& Qli_prev = Qlc.at(detail::mvaoi_key(S, Nr));
            for (std::size_t j = 0; j < J; ++j) {
                T s = zero;
                for (std::size_t t = 0; t < R; ++t) s += Qli_prev(j, t);
                betaLI(j, r) = Dli(j, r) * (one + s);
            }
        }

        // population conservation A X = Nn over the classes with Nn_r > 0
        std::vector<std::size_t> idx;
        for (std::size_t r = 0; r < R; ++r)
            if (Nn[r] > 0) idx.push_back(r);
        const std::size_t m = idx.size();
        Matrix<T> A(m, m, zero);
        for (std::size_t a = 0; a < m; ++a) {
            const std::size_t r = idx[a];
            for (std::size_t b = 0; b < m; ++b) {
                const std::size_t s = idx[b];
                T val = zero;
                if (s == r) {
                    val = Z[r];
                    for (std::size_t j = 0; j < J; ++j) val += betaLI(j, r);
                    for (std::size_t i = 0; i < K; ++i) val += Dt(i, r) * (one + Qsub[i](r, r));
                } else {
                    for (std::size_t i = 0; i < K; ++i) val += Dt(i, s) * Qsub[i](s, r);
                }
                A(a, b) = val;
            }
        }

        std::vector<T> Xk(R, zero);
        for (std::size_t a = 0; a < m; ++a) {
            const std::size_t r = idx[a];
            T denom = A(a, a);
            std::vector<int> Nr = Nn;
            Nr[r] -= 1;
            const std::vector<T>& Xner = Xc.at(detail::mvaoi_key(S, Nr));
            for (std::size_t b = 0; b < m; ++b) {
                if (b == a) continue;
                const std::size_t s = idx[b];
                std::vector<int> Ns = Nn;
                Ns[s] -= 1;
                const std::vector<T>& Xnes = Xc.at(detail::mvaoi_key(S, Ns));
                if (Xnes[r] > zero) denom += A(a, b) * (Xner[s] / Xnes[r]);
            }
            if (denom > zero) Xk[r] = num_traits<T>::from_int(Nn[r]) / denom;
        }

        Matrix<T> Qk_li(J, R, zero);
        for (std::size_t r = 0; r < R; ++r) {
            if (Nn[r] == 0) continue;
            for (std::size_t j = 0; j < J; ++j) Qk_li(j, r) = Xk[r] * betaLI(j, r);
        }
        for (std::size_t i = 0; i < K; ++i) {
            std::vector<T> U(R, zero), Qi(R, zero), Drow(R, zero);
            for (std::size_t r = 0; r < R; ++r) {
                Drow[r] = Dt(i, r);
                U[r] = Dt(i, r) * Xk[r];
            }
            for (std::size_t r = 0; r < R; ++r) {
                T s = U[r];
                for (std::size_t t = 0; t < R; ++t) s += U[t] * Qsub[i](t, r);
                Qi[r] = s;
            }
            Qc[i][key] = Qi;
            Dc[i][key] = Drow;
        }
        Xc[key] = Xk;
        Qlc[key] = Qk_li;
    }

    const std::vector<int> keyN = detail::mvaoi_key(zeroS, N);
    MvaoiResult<T> res;
    res.X = Xc.at(keyN);
    res.Qoi = Matrix<T>(K, R, zero);
    for (std::size_t i = 0; i < K; ++i) {
        const std::vector<T>& q = Qc[i].at(keyN);
        for (std::size_t r = 0; r < R; ++r) res.Qoi(i, r) = q[r];
    }
    res.Qli = Qlc.at(keyN);
    res.Qdelay.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r) res.Qdelay[r] = res.X[r] * Z[r];
    if (!want_soi) return res;

    // ---- mean number of in-service jobs per class at each OI station ----------------
    std::vector<std::size_t> shp(R), stride(R, 1);
    std::size_t total = 1;
    for (std::size_t d = 0; d < R; ++d) shp[d] = static_cast<std::size_t>(N[d]) + 1;
    for (std::size_t d = 1; d < R; ++d) stride[d] = stride[d - 1] * shp[d - 1];
    for (std::size_t d = 0; d < R; ++d) total *= shp[d];
    std::vector<std::vector<int>> subs(total, std::vector<int>(R, 0));
    std::vector<int> ssum(total, 0);
    for (std::size_t i = 0; i < total; ++i) {
        std::size_t li = i;
        for (std::size_t d = 0; d < R; ++d) {
            subs[i][d] = static_cast<int>(li % shp[d]);
            li /= shp[d];
            ssum[i] += subs[i][d];
        }
    }
    std::vector<std::size_t> ord(total);
    for (std::size_t i = 0; i < total; ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(),
                     [&](std::size_t a, std::size_t b) { return ssum[a] < ssum[b]; });

    Matrix<T> Xlat(total, R, zero);
    for (std::size_t i = 0; i < total; ++i) {
        const std::vector<T>& x = Xc.at(detail::mvaoi_key(zeroS, subs[i]));
        for (std::size_t r = 0; r < R; ++r) Xlat(i, r) = x[r];
    }

    res.Soi = Matrix<T>(K, R, zero);
    for (std::size_t mm = 0; mm < K; ++mm) {
        const OiInsvcResult<T> gm = pfqn_oi_insvc<T>(mu[mm], N);
        std::vector<T> muv(total, zero);
        for (std::size_t i = 0; i < total; ++i)
            if (ssum[i] > 0) muv[i] = mu[mm](subs[i]);
        Matrix<T> pMv(total, total, zero);
        pMv(0, 0) = one;
        for (std::size_t bb = 0; bb < total; ++bb) {
            const std::size_t b = ord[bb];
            if (ssum[b] == 0) continue;
            T acc0 = zero;
            for (std::size_t aa = 0; aa < total; ++aa) {
                const std::size_t a = ord[aa];
                if (ssum[a] == 0) continue;
                bool fits = true;
                for (std::size_t r = 0; r < R; ++r)
                    if (subs[a][r] > subs[b][r]) {
                        fits = false;
                        break;
                    }
                if (!fits || !(muv[a] > zero)) continue;
                T acc = zero;
                for (std::size_t r = 0; r < R; ++r)
                    if (subs[a][r] > 0)
                        acc += Xlat(b, r) * pMv(a - stride[r], b - stride[r]);
                pMv(a, b) = acc / muv[a];
                acc0 += pMv(a, b);
            }
            pMv(0, b) = one - acc0;  // empty-state probability by complement
        }
        const std::size_t idxN = total - 1;
        for (std::size_t r = 0; r < R; ++r) {
            T s = zero;
            for (std::size_t a = 0; a < total; ++a) s += pMv(a, idxN) * gm.g(a, r);
            res.Soi(mm, r) = s;
        }
    }
    return res;
}

/** Overload with unit visits. */
template <class T>
MvaoiResult<T> pfqn_mvaoi(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli, bool want_soi) {
    return pfqn_mvaoi(Z, N, mu, Dli, Matrix<T>(), want_soi);
}

/** Overload without the in-service means. */
template <class T>
MvaoiResult<T> pfqn_mvaoi(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli, const Matrix<T>& visits) {
    return pfqn_mvaoi(Z, N, mu, Dli, visits, false);
}

/** Overload without the in-service means, unit visits. */
template <class T>
MvaoiResult<T> pfqn_mvaoi(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli) {
    return pfqn_mvaoi(Z, N, mu, Dli, Matrix<T>(), false);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MVAOI_H
