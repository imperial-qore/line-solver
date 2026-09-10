/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_DESCRIPTOR_H
#define LINE_API_MDD_MDD_DESCRIPTOR_H

/**
 * Kronecker rate descriptor of a single-class closed queueing network.
 *
 * Port of matlab/src/api/mdd/mdd_descriptor.m, jline.api.mdd.Mdd_descriptor and
 * python/line_solver/api/mdd/descriptor.py, for the Miner-Ciardo-Donatelli
 * aggregation `mdd_mcd` (SIGMETRICS 2000).
 *
 * The transition rate matrix is expressed compositionally as
 * R = sum_e (kron_k W_k^e) restricted to the reachable set, with
 * W_k^e[i,j] = lambda_k^e[i] * Prob_k^e(i,j) (Eq. 1). Each level k is a station
 * and each event e is a completion at station a routed to b.
 *
 * EXPONENTIAL STATIONS. The local state is the population alone:
 * W_a^e[i,i-1] = mu[a]*min(i,servers[a])*P[a][b] for i >= 1 (departure),
 * W_b^e[i,i+1] = 1 for i <= N-1 (arrival), and the identity elsewhere.
 *
 * PHASE-TYPE STATIONS. The local state is the PAIR (population, phase of the
 * job in service), encoded in one level rather than two. Splitting them does not
 * work: on a completion routed into station b the phase at b restarts only when
 * b was empty, a joint condition on b's two components, which is not a product
 * of per-level terms. Merging them keeps every event local:
 *
 *   index 0                   : station empty
 *   index 1 + (n-1)*h + (a-1) : n jobs present, job in service in phase a
 *   domain                    = 1 + N*h    (h = 1 reproduces index = n)
 *
 * With exit vector t = D1*1 and entry law pie, departure (n,a)->(n-1,b) at
 * t[a]*P[a][b]*pie[b] for n >= 2 and (1,a)->0 at t[a]*P[a][b]; arrival 0->(1,b)
 * at pie[b] and (m,c)->(m+1,c) at 1 for m >= 1; internal (n,a)->(n,b) at
 * D0[a][b] for n >= 1, a != b.
 *
 * RESTRICTIONS. A phase-type station must be single-server: with c > 1 or an
 * infinite server the local state would have to count jobs per phase rather than
 * name one phase, a different and much larger encoding. It must also be
 * NON-preemptive, because the composite level names the phase of the one job in
 * service and restarts it at pie when the next job starts; under preemptive
 * resume an arrival suspends that job and its phase has to be remembered, so the
 * local state would need a stack of phases. That matters for LCFSPR, which is
 * BCMP type 2 and stays product-form under general service: the insensitivity is
 * real but is NOT reachable through this encoding. Exponential service is
 * unaffected, preemption being immaterial by memorylessness. Pass the
 * disciplines to have the case rejected rather than silently modelled as
 * non-preemptive.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mdd/mdd_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mdd {

namespace detail {

/** Disciplines whose preemptive-resume semantics the composite encoding cannot carry. */
inline bool desc_is_preemptive_resume(const std::string& nm) {
    return nm == "LCFSPR" || nm == "FCFSPR" || nm == "LCFSPRPRIO" || nm == "FCFSPRPRIO";
}

/** Shared-server disciplines, where every job present holds its own phase. */
inline bool desc_is_shared_server(const std::string& nm) {
    return nm == "PS" || nm == "DPS" || nm == "GPS";
}

inline std::string desc_upper(const std::string& s) {
    std::string out = s;
    for (std::size_t i = 0; i < out.size(); ++i)
        if (out[i] >= 'a' && out[i] <= 'z') out[i] = static_cast<char>(out[i] - 'a' + 'A');
    return out;
}

/** Local index of (population n, service phase a); 0 when the station is empty. */
inline int desc_idx(int n, int a, std::size_t h) {
    if (n == 0) return 0;
    return 1 + (n - 1) * static_cast<int>(h) + (a - 1);
}

/** Population of a local index; 0 when empty. */
inline int desc_population(int index, std::size_t h) {
    if (index == 0) return 0;
    return (index - 1) / static_cast<int>(h) + 1;
}

/** Service phase (1-based) of a local index; 0 when empty. */
inline int desc_phase(int index, std::size_t h) {
    if (index == 0) return 0;
    return (index - 1) % static_cast<int>(h) + 1;
}

/** Phase changes that do not complete a service, at any population n >= 1. */
template <class T>
MddLocalMatrix<T> desc_internal(const Matrix<T>& D0i, std::size_t h, int N, int d) {
    const T zero = num_traits<T>::from_int(0);
    typename MddLocalMatrix<T>::Builder bld(static_cast<std::size_t>(d));
    for (int n = 1; n <= N; ++n)
        for (std::size_t a = 1; a <= h; ++a)
            for (std::size_t b = 1; b <= h; ++b) {
                if (a == b || D0i(a - 1, b - 1) == zero) continue;
                bld.add(static_cast<std::size_t>(desc_idx(n, static_cast<int>(a), h)),
                        static_cast<std::size_t>(desc_idx(n, static_cast<int>(b), h)),
                        D0i(a - 1, b - 1));
            }
    return bld.build();
}

/** Completion at this station, routed out with probability pr. */
template <class T>
MddLocalMatrix<T> desc_departure(const Matrix<T>& D1i, const std::vector<T>& piei, std::size_t h,
                                 int N, int d, const T& mui, double srv, const T& pr) {
    const T zero = num_traits<T>::from_int(0);
    typename MddLocalMatrix<T>::Builder bld(static_cast<std::size_t>(d));
    if (h == 1) {
        // exponential: the multi-server and delay rate laws live here
        for (int n = 1; n <= N; ++n) {
            const double cap = static_cast<double>(n) < srv ? static_cast<double>(n) : srv;
            bld.add(static_cast<std::size_t>(n), static_cast<std::size_t>(n - 1),
                    T(mui * num_traits<T>::from_double(cap) * pr));
        }
    } else {
        std::vector<T> t(h, zero);
        for (std::size_t a = 0; a < h; ++a) {
            T s = zero;
            for (std::size_t b = 0; b < h; ++b) s += D1i(a, b);
            t[a] = s;
        }
        for (int n = 1; n <= N; ++n)
            for (std::size_t a = 1; a <= h; ++a) {
                if (t[a - 1] == zero) continue;
                if (n == 1) {
                    bld.add(static_cast<std::size_t>(desc_idx(1, static_cast<int>(a), h)), 0,
                            T(t[a - 1] * pr));
                } else {
                    for (std::size_t b = 1; b <= h; ++b) {
                        if (piei[b - 1] == zero) continue;
                        bld.add(static_cast<std::size_t>(desc_idx(n, static_cast<int>(a), h)),
                                static_cast<std::size_t>(
                                    desc_idx(n - 1, static_cast<int>(b), h)),
                                T(t[a - 1] * pr * piei[b - 1]));
                    }
                }
            }
    }
    return bld.build();
}

/** An arrival starts service only when the station was empty. */
template <class T>
MddLocalMatrix<T> desc_arrival(const std::vector<T>& pieb, std::size_t h, int N, int d) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    typename MddLocalMatrix<T>::Builder bld(static_cast<std::size_t>(d));
    if (h == 1) {
        for (int m = 0; m < N; ++m)
            bld.add(static_cast<std::size_t>(m), static_cast<std::size_t>(m + 1), one);
    } else {
        for (std::size_t b = 1; b <= h; ++b) {
            if (pieb[b - 1] == zero) continue;
            bld.add(0, static_cast<std::size_t>(desc_idx(1, static_cast<int>(b), h)), pieb[b - 1]);
        }
        for (int m = 1; m < N; ++m)
            for (std::size_t c = 1; c <= h; ++c)
                bld.add(static_cast<std::size_t>(desc_idx(m, static_cast<int>(c), h)),
                        static_cast<std::size_t>(desc_idx(m + 1, static_cast<int>(c), h)), one);
    }
    return bld.build();
}

/** Successor local-index vectors of a state, used to generate the reachable set. */
template <class T>
std::vector<std::vector<int>> desc_successors(const std::vector<int>& s,
                                              const std::vector<Matrix<T>>& D0,
                                              const std::vector<Matrix<T>>& D1,
                                              const std::vector<std::vector<T>>& pie,
                                              const std::vector<std::size_t>& h,
                                              const Matrix<T>& P) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = s.size();
    std::vector<std::vector<int>> out;
    for (std::size_t i = 0; i < K; ++i) {
        const int ni = desc_population(s[i], h[i]);
        const int ai = desc_phase(s[i], h[i]);
        if (ni == 0) continue;
        // internal phase change
        if (h[i] > 1) {
            for (std::size_t b = 1; b <= h[i]; ++b) {
                if (static_cast<int>(b) == ai || D0[i](ai - 1, b - 1) == zero) continue;
                std::vector<int> t = s;
                t[i] = desc_idx(ni, static_cast<int>(b), h[i]);
                out.push_back(t);
            }
        }
        // completion routed to j
        T exits = one;
        if (h[i] > 1) {
            exits = zero;
            for (std::size_t b = 0; b < h[i]; ++b) exits += D1[i](ai - 1, b);
        }
        if (exits == zero) continue;
        for (std::size_t j = 0; j < K; ++j) {
            if (j == i || !(P(i, j) > zero)) continue;
            const int nj = desc_population(s[j], h[j]);
            const int aj = desc_phase(s[j], h[j]);
            int newi = (h[i] == 1) ? desc_idx(ni - 1, 1, 1) : 0;
            for (std::size_t bi = 1; bi <= h[i]; ++bi) {
                if (h[i] > 1) {
                    if (ni == 1) {
                        newi = 0;
                    } else if (pie[i][bi - 1] == zero) {
                        continue;
                    } else {
                        newi = desc_idx(ni - 1, static_cast<int>(bi), h[i]);
                    }
                } else if (bi > 1) {
                    continue;
                }
                for (std::size_t bj = 1; bj <= h[j]; ++bj) {
                    int newj;
                    if (nj == 0) {
                        if (pie[j][bj - 1] == zero) continue;
                        newj = desc_idx(1, static_cast<int>(bj), h[j]);
                    } else if (bj > 1) {
                        continue;
                    } else {
                        newj = desc_idx(nj + 1, aj, h[j]);
                    }
                    std::vector<int> t = s;
                    t[i] = newi;
                    t[j] = newj;
                    out.push_back(t);
                }
                if (ni == 1 && h[i] > 1) break;
            }
        }
    }
    return out;
}

}  // namespace detail

/**
 * Build the descriptor.
 *
 * @param mu station service rates, 1/E[S]; entry i is ignored when station i is
 *        given a phase-type law through proc
 * @param P station-to-station routing matrix, row-stochastic
 * @param servers servers per station, infinite for delay/IS
 * @param N closed population
 * @param proc per-station service law; an absent or `present == false` entry is
 *        an exponential station
 * @param sched per-station discipline names, consulted only to REJECT a
 *        phase-type law at a preemptive-resume or shared-server station; may be
 *        empty when every station is non-preemptive
 */
template <class T>
MddDescriptor<T> mdd_descriptor(const std::vector<T>& mu, const Matrix<T>& P,
                                const std::vector<double>& servers, int N,
                                const std::vector<MddServiceLaw<T>>& proc =
                                    std::vector<MddServiceLaw<T>>(),
                                const std::vector<std::string>& sched =
                                    std::vector<std::string>()) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = mu.size();
    if (K == 0) throw InputError("mdd_descriptor: the network has no stations");
    if (P.rows() != K || P.cols() != K)
        throw InputError("mdd_descriptor: the routing matrix is not (K x K)");
    if (servers.size() != K)
        throw InputError("mdd_descriptor: one server count per station is required");
    if (N < 0) throw InputError("mdd_descriptor: the population must be non-negative");

    std::vector<Matrix<T>> D0(K), D1(K);
    std::vector<std::vector<T>> entry(K);
    std::vector<std::size_t> h(K, 1);

    for (std::size_t i = 0; i < K; ++i) {
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
        entry[i] = mdd_entry_law(proc[i].pie, D1[i], h[i], i, "mdd_descriptor");
        if (h[i] > 1 && servers[i] != 1)
            throw InputError("mdd_descriptor: station " + std::to_string(i + 1) +
                             " has a phase-type service law and " + std::to_string(servers[i]) +
                             " servers; a multi-server or delay station would have to count jobs "
                             "per phase rather than name the phase of one job in service, which "
                             "this encoding does not carry");
        if (h[i] > 1 && i < sched.size() && !sched[i].empty()) {
            const std::string nm = detail::desc_upper(sched[i]);
            if (detail::desc_is_preemptive_resume(nm))
                throw InputError("mdd_descriptor: station " + std::to_string(i + 1) +
                                 " combines a phase-type service law with a preemptive-resume "
                                 "discipline; the suspended jobs' phases would have to be stacked "
                                 "in the local state, which this encoding does not carry, and the "
                                 "descriptor would silently model the non-preemptive chain "
                                 "instead");
            if (detail::desc_is_shared_server(nm))
                throw InputError("mdd_descriptor: station " + std::to_string(i + 1) +
                                 " combines a phase-type service law with a shared-server "
                                 "discipline; every job present is in service and holds its own "
                                 "phase, which this encoding does not carry. Use mdd_ps, whose "
                                 "local state is the per-phase count vector.");
        }
    }

    std::vector<int> d(K, 0);
    for (std::size_t i = 0; i < K; ++i) d[i] = 1 + N * static_cast<int>(h[i]);

    MddDescriptor<T> desc;
    desc.K = K;
    desc.N = N;
    desc.domain = d;
    desc.mu = mu;
    desc.servers = servers;
    desc.P = P;
    desc.nphases = h;

    // ---- index maps
    desc.valuemap.assign(K, std::vector<double>());
    for (std::size_t i = 0; i < K; ++i) {
        desc.valuemap[i].assign(static_cast<std::size_t>(d[i]), 0.0);
        for (int n = 1; n <= N; ++n)
            for (std::size_t a = 0; a < h[i]; ++a)
                desc.valuemap[i][static_cast<std::size_t>(1 + (n - 1) * static_cast<int>(h[i]) +
                                                          static_cast<int>(a))] =
                    static_cast<double>(n);
    }

    // ---- initial state: all jobs at station 1, in its entry phase
    desc.init.assign(K, 0);
    std::size_t first_phase = 0;
    for (std::size_t a = 0; a < h[0]; ++a)
        if (entry[0][a] > zero) {
            first_phase = a;
            break;
        }
    desc.init[0] = detail::desc_idx(N, static_cast<int>(first_phase) + 1, h[0]);

    const std::vector<Matrix<T>> fD0 = D0, fD1 = D1;
    const std::vector<std::vector<T>> fpie = entry;
    const std::vector<std::size_t> fh = h;
    const Matrix<T> fP = P;
    desc.nextfun = [fD0, fD1, fpie, fh, fP](const std::vector<int>& state) {
        return detail::desc_successors(state, fD0, fD1, fpie, fh, fP);
    };

    // ---- events: internal phase changes, one per phase-type station
    for (std::size_t i = 0; i < K; ++i) {
        if (h[i] == 1) continue;
        const MddLocalMatrix<T> Wi = detail::desc_internal(D0[i], h[i], N, d[i]);
        if (Wi.nnz == 0) continue;
        MddEvent<T> ev;
        ev.a = i;
        ev.b = i;
        ev.lev.push_back(i);
        ev.W.push_back(Wi);
        desc.events.push_back(ev);
    }
    // ---- events: completions routed a -> b
    for (std::size_t a = 0; a < K; ++a)
        for (std::size_t b = 0; b < K; ++b) {
            if (a == b || !(P(a, b) > zero)) continue;
            MddEvent<T> ev;
            ev.a = a;
            ev.b = b;
            ev.lev.push_back(a);
            ev.lev.push_back(b);
            ev.W.push_back(detail::desc_departure(D1[a], entry[a], h[a], N, d[a], mu[a],
                                                  servers[a], P(a, b)));
            ev.W.push_back(detail::desc_arrival(entry[b], h[b], N, d[b]));
            desc.events.push_back(ev);
        }
    return desc;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_DESCRIPTOR_H
