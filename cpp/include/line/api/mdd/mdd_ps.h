/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_PS_H
#define LINE_API_MDD_MDD_PS_H

/**
 * Kronecker rate descriptor for shared-server stations with phase-type service.
 *
 * Port of matlab/src/api/mdd/mdd_ps.m, jline.api.mdd.Mdd_ps and
 * python/line_solver/api/mdd/ps.py.
 *
 * Under processor sharing every job at a station is in service at once, each
 * holding its own phase, so naming a single in-service phase (what
 * `mdd_descriptor` does, which is non-preemptive semantics) cannot represent the
 * state. The local state here is instead the PER-PHASE COUNT vector
 * v = (v_1,...,v_h), v_a jobs in phase a, with n = sum(v) jobs present. That is
 * still a per-station quantity, so every event stays a product of per-level
 * terms and the Kronecker form of Eq. 1 survives.
 *
 * With one server shared by n jobs each job advances at rate 1/n, so from local
 * state v with n = sum(v):
 *
 *   internal   v -> v - e_a + e_b   at v_a * D0[a][b] / n     (a != b)
 *   departure  v -> v - e_a         at v_a * t[a] / n * P[i][j]
 *   arrival    v -> v + e_b         at pie[b]
 *
 * An infinite-server (delay) station is the same without the 1/n scaling. For
 * h = 1 the departure rate collapses to n*mu/n = mu at PS and to n*mu at IS,
 * reproducing the usual single-server and delay rate laws.
 *
 * The local domain is the number of compositions of 0..N over h phases,
 * C(N+h,h), against 1+N*h for the non-preemptive encoding: the price of tracking
 * every job's phase rather than one.
 */

#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "line/api/mdd/mdd_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mdd {

namespace detail {

/** Compositions of 0..N over h phases, in the reference's row order. */
inline std::vector<std::vector<int>> ps_compositions(std::size_t h, int N) {
    std::vector<std::vector<int>> out;
    if (h == 1) {
        for (int n = 0; n <= N; ++n) out.push_back(std::vector<int>(1, n));
        return out;
    }
    const std::vector<std::vector<int>> sub = ps_compositions(h - 1, N);
    for (int v1 = 0; v1 <= N; ++v1) {
        for (std::size_t r = 0; r < sub.size(); ++r) {
            int s = 0;
            for (std::size_t a = 0; a < sub[r].size(); ++a) s += sub[r][a];
            if (s > N - v1) continue;
            std::vector<int> row;
            row.push_back(v1);
            row.insert(row.end(), sub[r].begin(), sub[r].end());
            out.push_back(row);
        }
    }
    return out;
}

/** Service-rate scaling: 1/n shared by n jobs at PS, unscaled at IS. */
template <class T>
T ps_share(int n, double srv) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (n == 0) return zero;
    if (std::isinf(srv)) return one;
    return T(one / num_traits<T>::from_int(n));
}

inline std::size_t ps_lookup(const std::map<std::vector<int>, std::size_t>& lut,
                             const std::vector<int>& v) {
    const std::map<std::vector<int>, std::size_t>::const_iterator it = lut.find(v);
    if (it == lut.end())
        throw InputError("mdd_ps: a per-phase count vector left the composition state space");
    return it->second;
}

template <class T>
MddLocalMatrix<T> ps_internal(const Matrix<T>& D0i, const std::vector<std::vector<int>>& C,
                              const std::map<std::vector<int>, std::size_t>& lut, std::size_t h,
                              std::size_t d, double srv) {
    const T zero = num_traits<T>::from_int(0);
    typename MddLocalMatrix<T>::Builder bld(d);
    for (std::size_t r = 0; r < d; ++r) {
        const std::vector<int>& v = C[r];
        int n = 0;
        for (std::size_t a = 0; a < h; ++a) n += v[a];
        if (n == 0) continue;
        const T sc = ps_share<T>(n, srv);
        for (std::size_t a = 0; a < h; ++a) {
            if (v[a] == 0) continue;
            for (std::size_t b = 0; b < h; ++b) {
                if (a == b || D0i(a, b) == zero) continue;
                std::vector<int> w = v;
                --w[a];
                ++w[b];
                bld.add(r, ps_lookup(lut, w),
                        T(num_traits<T>::from_int(v[a]) * D0i(a, b) * sc));
            }
        }
    }
    return bld.build();
}

template <class T>
MddLocalMatrix<T> ps_departure(const Matrix<T>& D1i, const std::vector<std::vector<int>>& C,
                               const std::map<std::vector<int>, std::size_t>& lut, std::size_t h,
                               std::size_t d, double srv, const T& pr) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> t(h, zero);
    for (std::size_t a = 0; a < h; ++a) {
        T s = zero;
        for (std::size_t b = 0; b < h; ++b) s += D1i(a, b);
        t[a] = s;
    }
    typename MddLocalMatrix<T>::Builder bld(d);
    for (std::size_t r = 0; r < d; ++r) {
        const std::vector<int>& v = C[r];
        int n = 0;
        for (std::size_t a = 0; a < h; ++a) n += v[a];
        if (n == 0) continue;
        const T sc = ps_share<T>(n, srv);
        for (std::size_t a = 0; a < h; ++a) {
            if (v[a] == 0 || t[a] == zero) continue;
            std::vector<int> w = v;
            --w[a];
            bld.add(r, ps_lookup(lut, w),
                    T(num_traits<T>::from_int(v[a]) * t[a] * sc * pr));
        }
    }
    return bld.build();
}

template <class T>
MddLocalMatrix<T> ps_arrival(const std::vector<T>& pieb, const std::vector<std::vector<int>>& C,
                             const std::map<std::vector<int>, std::size_t>& lut, std::size_t h,
                             int N, std::size_t d) {
    const T zero = num_traits<T>::from_int(0);
    typename MddLocalMatrix<T>::Builder bld(d);
    for (std::size_t r = 0; r < d; ++r) {
        const std::vector<int>& v = C[r];
        int n = 0;
        for (std::size_t a = 0; a < h; ++a) n += v[a];
        if (n >= N) continue;
        for (std::size_t b = 0; b < h; ++b) {
            if (pieb[b] == zero) continue;
            std::vector<int> w = v;
            ++w[b];
            bld.add(r, ps_lookup(lut, w), pieb[b]);
        }
    }
    return bld.build();
}

template <class T>
std::vector<std::vector<int>> ps_successors(
    const std::vector<int>& s, const std::vector<Matrix<T>>& D0, const std::vector<Matrix<T>>& D1,
    const std::vector<std::vector<T>>& pie, const std::vector<std::size_t>& h,
    const std::vector<std::vector<std::vector<int>>>& comp,
    const std::vector<std::map<std::vector<int>, std::size_t>>& lut, int N, const Matrix<T>& P) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = s.size();
    std::vector<std::vector<int>> out;
    for (std::size_t i = 0; i < K; ++i) {
        const std::vector<int>& v = comp[i][static_cast<std::size_t>(s[i])];
        int n = 0;
        for (std::size_t a = 0; a < h[i]; ++a) n += v[a];
        if (n == 0) continue;
        // internal phase moves
        for (std::size_t a = 0; a < h[i]; ++a) {
            if (v[a] == 0) continue;
            for (std::size_t b = 0; b < h[i]; ++b) {
                if (a == b || D0[i](a, b) == zero) continue;
                std::vector<int> w = v;
                --w[a];
                ++w[b];
                std::vector<int> t = s;
                t[i] = static_cast<int>(ps_lookup(lut[i], w));
                out.push_back(t);
            }
        }
        // completions routed to j
        for (std::size_t a = 0; a < h[i]; ++a) {
            if (v[a] == 0) continue;
            T ta = zero;
            for (std::size_t b = 0; b < h[i]; ++b) ta += D1[i](a, b);
            if (ta == zero) continue;
            std::vector<int> w = v;
            --w[a];
            for (std::size_t j = 0; j < K; ++j) {
                if (j == i || !(P(i, j) > zero)) continue;
                const std::vector<int>& vj = comp[j][static_cast<std::size_t>(s[j])];
                int nj = 0;
                for (std::size_t b = 0; b < h[j]; ++b) nj += vj[b];
                if (nj >= N) continue;
                for (std::size_t b = 0; b < h[j]; ++b) {
                    if (pie[j][b] == zero) continue;
                    std::vector<int> wj = vj;
                    ++wj[b];
                    std::vector<int> t = s;
                    t[i] = static_cast<int>(ps_lookup(lut[i], w));
                    t[j] = static_cast<int>(ps_lookup(lut[j], wj));
                    out.push_back(t);
                }
            }
        }
    }
    return out;
}

}  // namespace detail

/**
 * Build the descriptor.
 *
 * @param mu station service rates, ignored where proc gives a law
 * @param P station-to-station routing matrix, row-stochastic
 * @param servers servers per station, 1 (PS) or infinite (IS); no other value
 *        has a per-phase-count encoding here
 * @param N closed population
 * @param proc per-station service law; an absent or `present == false` entry is
 *        an exponential station
 */
template <class T>
MddDescriptor<T> mdd_ps(const std::vector<T>& mu, const Matrix<T>& P,
                        const std::vector<double>& servers, int N,
                        const std::vector<MddServiceLaw<T>>& proc =
                            std::vector<MddServiceLaw<T>>()) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = mu.size();
    if (K == 0) throw InputError("mdd_ps: the network has no stations");
    if (P.rows() != K || P.cols() != K)
        throw InputError("mdd_ps: the routing matrix is not (K x K)");
    if (servers.size() != K) throw InputError("mdd_ps: one server count per station is required");

    std::vector<Matrix<T>> D0(K), D1(K);
    std::vector<std::vector<T>> entry(K);
    std::vector<std::size_t> h(K, 1);

    for (std::size_t i = 0; i < K; ++i) {
        if (!(servers[i] == 1 || std::isinf(servers[i])))
            throw InputError("mdd_ps: station " + std::to_string(i + 1) + " has " +
                             std::to_string(servers[i]) +
                             " servers; only processor sharing (1) and infinite server have a "
                             "per-phase-count encoding here");
        const bool has_law = i < proc.size() && proc[i].present;
        if (!has_law) {
            D0[i] = Matrix<T>(1, 1, T(-mu[i]));
            D1[i] = Matrix<T>(1, 1, mu[i]);
            entry[i] = std::vector<T>(1, one);
            h[i] = 1;
            continue;
        }
        D0[i] = proc[i].D0;
        D1[i] = proc[i].D1;
        h[i] = proc[i].phases();
        entry[i] = mdd_entry_law(proc[i].pie, D1[i], h[i], i, "mdd_ps");
    }

    // ---- per-station composition state space
    std::vector<std::vector<std::vector<int>>> comp(K);
    std::vector<std::map<std::vector<int>, std::size_t>> lut(K);
    std::vector<int> d(K, 0);
    for (std::size_t i = 0; i < K; ++i) {
        comp[i] = detail::ps_compositions(h[i], N);
        for (std::size_t r = 0; r < comp[i].size(); ++r) lut[i][comp[i][r]] = r;
        d[i] = static_cast<int>(comp[i].size());
    }

    MddDescriptor<T> desc;
    desc.K = K;
    desc.N = N;
    desc.domain = d;
    desc.mu = mu;
    desc.servers = servers;
    desc.P = P;
    desc.nphases = h;
    desc.valuemap.assign(K, std::vector<double>());
    for (std::size_t i = 0; i < K; ++i) {
        desc.valuemap[i].assign(static_cast<std::size_t>(d[i]), 0.0);
        for (std::size_t r = 0; r < static_cast<std::size_t>(d[i]); ++r) {
            int s = 0;
            for (std::size_t a = 0; a < h[i]; ++a) s += comp[i][r][a];
            desc.valuemap[i][r] = static_cast<double>(s);
        }
    }

    // ---- initial state: all jobs at station 1, entered in its entry phase
    desc.init.assign(K, 0);
    for (std::size_t i = 0; i < K; ++i) {
        std::vector<int> v(h[i], 0);
        if (i == 0) {
            std::size_t first = 0;
            for (std::size_t a = 0; a < h[0]; ++a)
                if (entry[0][a] > zero) {
                    first = a;
                    break;
                }
            v[first] = N;
        }
        desc.init[i] = static_cast<int>(detail::ps_lookup(lut[i], v));
    }

    const std::vector<Matrix<T>> fD0 = D0, fD1 = D1;
    const std::vector<std::vector<T>> fpie = entry;
    const std::vector<std::size_t> fh = h;
    const std::vector<std::vector<std::vector<int>>> fcomp = comp;
    const std::vector<std::map<std::vector<int>, std::size_t>> flut = lut;
    const Matrix<T> fP = P;
    const int fN = N;
    desc.nextfun = [fD0, fD1, fpie, fh, fcomp, flut, fN, fP](const std::vector<int>& state) {
        return detail::ps_successors(state, fD0, fD1, fpie, fh, fcomp, flut, fN, fP);
    };

    // ---- events
    for (std::size_t i = 0; i < K; ++i) {
        if (h[i] == 1) continue;
        const MddLocalMatrix<T> Wi = detail::ps_internal(D0[i], comp[i], lut[i], h[i],
                                                         static_cast<std::size_t>(d[i]),
                                                         servers[i]);
        if (Wi.nnz == 0) continue;
        MddEvent<T> ev;
        ev.a = i;
        ev.b = i;
        ev.lev.push_back(i);
        ev.W.push_back(Wi);
        desc.events.push_back(ev);
    }
    for (std::size_t a = 0; a < K; ++a)
        for (std::size_t b = 0; b < K; ++b) {
            if (a == b || !(P(a, b) > zero)) continue;
            MddEvent<T> ev;
            ev.a = a;
            ev.b = b;
            ev.lev.push_back(a);
            ev.lev.push_back(b);
            ev.W.push_back(detail::ps_departure(D1[a], comp[a], lut[a], h[a],
                                                static_cast<std::size_t>(d[a]), servers[a],
                                                P(a, b)));
            ev.W.push_back(detail::ps_arrival(entry[b], comp[b], lut[b], h[b], N,
                                              static_cast<std::size_t>(d[b])));
            desc.events.push_back(ev);
        }
    return desc;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_PS_H
