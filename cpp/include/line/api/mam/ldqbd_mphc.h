/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_LDQBD_MPHC_H
#define LINE_API_MAM_LDQBD_MPHC_H

/**
 * @file
 * @ingroup api_mam
 * Port of `ldqbd_mphc.m` and `ph_multisets.m`: the exact level-dependent QBD
 * blocks of an M/PH/c queue.
 *
 * THE POINT OF THE MULTISET. The level is the number of jobs at the station,
 * and the coordinate INSIDE a level is the multiset of the phases the min(n,c)
 * busy servers sit in. The collapsed alternative -- one PH process run at
 * min(n,c) times its speed -- gets the aggregate service rate right but forgets
 * which phase each busy server is in, which turns c servers into one fast
 * server whose remaining work is a single phase-type variable. Counting rather
 * than ordering the phases costs nchoosek(min(n,c)+p-1, p-1) states per level
 * instead of p^min(n,c), because identical servers are exchangeable.
 *
 * LEVEL SIZES GROW over the boundary levels 0..c and repeat above them, so the
 * blocks joining differently sized neighbours are rectangular. `ldqbd`, its
 * rate matrices and its stationary vector all accept that heterogeneity; level
 * 0 is the single empty configuration.
 *
 * References: S. Asmussen and J.R. Moller, "Calculation of the steady state
 * waiting time distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems
 * 37(1):9-29, 2001; M. F. Neuts, "Matrix-geometric solutions in stochastic
 * models", Johns Hopkins University Press, 1981.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * The widest level is the repeating one, and the LD-QBD recursion inverts one
 * matrix of that order per level, so that is the size worth guarding.
 */
inline constexpr std::size_t LDQBD_MPHC_MAX_CONFIGS = 2000;

/**
 * Configurations of k identical servers over p service phases.
 *
 * Rows are the compositions of k into p nonnegative parts: entry (r,i) is the
 * number of the k busy servers sitting in phase i. The order is fixed and
 * shared by every caller, so a configuration index means the same thing in each
 * of them: the first part descends, hence k = 1 yields the identity rows
 * e_1 ... e_p in phase order, which is what makes the c = 1 case coincide with
 * plain phase indexing.
 */
inline std::vector<std::vector<int> > ph_multisets(std::size_t p, std::size_t k) {
    std::vector<std::vector<int> > out;
    if (k == 0) {
        out.push_back(std::vector<int>(p, 0));
        return out;
    }
    if (p == 1) {
        out.push_back(std::vector<int>(1, static_cast<int>(k)));
        return out;
    }
    for (int first = static_cast<int>(k); first >= 0; --first) {
        const std::vector<std::vector<int> > rest =
            ph_multisets(p - 1, k - static_cast<std::size_t>(first));
        for (std::size_t i = 0; i < rest.size(); ++i) {
            std::vector<int> row;
            row.push_back(first);
            row.insert(row.end(), rest[i].begin(), rest[i].end());
            out.push_back(row);
        }
    }
    return out;
}

/** The three block lists of a level-dependent QBD, as `ldqbd` takes them. */
template <class T>
struct LdqbdMphcBlocks {
    std::vector<Matrix<T> > Q0;  // size Nlev   : upward, level n -> n+1
    std::vector<Matrix<T> > Q1;  // size Nlev+1 : local
    std::vector<Matrix<T> > Q2;  // size Nlev+1 : downward, Q2[n] leaves level n (Q2[0] unused)
};

/**
 * Carried-phase MAP/c blocks, the `carry` form of `ldqbd_mphc`. The station draws service
 * times from ONE MAP in start order: a job's service begins in phase h, the phase the previous
 * job's service ended in. V = (-D0)^-1 D1: a job starting in h ends in j w.p. V(h,j), so each
 * busy server is a pair (i,j) run under the Doob transform of D0 conditioned on ending in j
 * (i -> k at D0(i,k)V(k,j)/V(i,j), completion at D1(i,j)/V(i,j)). A level is the multiset of
 * busy pairs times h, the end phase of the job started last (index row*p + h); level 0 is h
 * alone. At c = 1 this is the frozen-phase MAP/MAP/1 chain; a renewal MAP gives the PH/c law.
 */
template <class T>
LdqbdMphcBlocks<T> ldqbd_mapc_carry(const Matrix<T>& D0, const Matrix<T>& D1, std::size_t cmax,
                                    const std::vector<T>& arrRate, const std::vector<T>& sf) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t p = D0.rows();
    const std::size_t Nlev = arrRate.size() - 1;
    Matrix<T> negD0(p, p, zero);
    for (std::size_t i = 0; i < p; ++i)
        for (std::size_t j = 0; j < p; ++j) negD0(i, j) = T(-D0(i, j));
    const Matrix<T> Ninv = inverse(negD0);
    Matrix<T> V(p, p, zero);
    std::vector<std::vector<int> > tp(p, std::vector<int>(p, -1));
    std::vector<std::pair<std::size_t, std::size_t> > pairs;
    for (std::size_t i = 0; i < p; ++i)
        for (std::size_t j = 0; j < p; ++j) {
            T v = zero;
            for (std::size_t k = 0; k < p; ++k) v += T(Ninv(i, k) * D1(k, j));
            if (std::fabs(num_traits<T>::to_double(v)) < 1e-14) v = zero;
            V(i, j) = v;
            if (num_traits<T>::to_double(v) > 0.0) {
                tp[i][j] = static_cast<int>(pairs.size());
                pairs.push_back(std::make_pair(i, j));
            }
        }
    const std::size_t NT = pairs.size();
    std::vector<T> done(NT, zero);
    for (std::size_t t = 0; t < NT; ++t)
        done[t] = T(D1(pairs[t].first, pairs[t].second) / V(pairs[t].first, pairs[t].second));

    std::vector<std::vector<std::vector<int> > > cfg(cmax + 1);
    std::vector<std::map<std::vector<int>, std::size_t> > pos(cmax + 1);
    std::vector<std::size_t> nCfg(cmax + 1, 0), nk(cmax + 1, 0);
    for (std::size_t k = 0; k <= cmax; ++k) {
        cfg[k] = ph_multisets(NT, k);
        for (std::size_t r = 0; r < cfg[k].size(); ++r) pos[k][cfg[k][r]] = r;
        nCfg[k] = cfg[k].size();
        nk[k] = nCfg[k] * p;
    }
    if (nk[cmax] > LDQBD_MPHC_MAX_CONFIGS)
        throw UnsupportedError(
            "ldqbd_mphc: the exact MAP/c chain with a carried service phase needs " +
            std::to_string(nk[cmax]) + " configurations per level for " + std::to_string(cmax) +
            " servers and " + std::to_string(p) + " service phases, above the " +
            std::to_string(LDQBD_MPHC_MAX_CONFIGS) +
            " the level-by-level inverses can carry. Use fewer servers or SolverLDES on this model");

    std::vector<Matrix<T> > LOC(cmax + 1), UP(cmax + 1), DN(cmax + 1);
    for (std::size_t k = 0; k <= cmax; ++k) {
        LOC[k] = Matrix<T>(nk[k], nk[k], zero);
        if (k < cmax) UP[k] = Matrix<T>(nk[k], nk[k + 1], zero);
        if (k > 0) DN[k] = Matrix<T>(nk[k], nk[k - 1], zero);
        for (std::size_t row = 0; row < nCfg[k]; ++row) {
            const std::vector<int>& m = cfg[k][row];
            for (std::size_t h = 0; h < p; ++h) {
                const std::size_t a = row * p + h;
                for (std::size_t t = 0; t < NT; ++t) {
                    if (m[t] == 0) continue;
                    const T mult = num_traits<T>::from_int(m[t]);
                    const std::size_t i = pairs[t].first, j = pairs[t].second;
                    LOC[k](a, a) += T(mult * D0(i, i));
                    for (std::size_t kk = 0; kk < p; ++kk) {
                        if (kk == i || tp[kk][j] < 0 || num_traits<T>::to_double(D0(i, kk)) == 0.0)
                            continue;
                        std::vector<int> mm = m;
                        --mm[t];
                        ++mm[tp[kk][j]];
                        LOC[k](a, pos[k][mm] * p + h) += T(mult * D0(i, kk) * V(kk, j) / V(i, j));
                    }
                    if (k > 0) {
                        std::vector<int> mm = m;
                        --mm[t];
                        DN[k](a, pos[k - 1][mm] * p + h) += T(mult * done[t]);
                    }
                }
                if (k < cmax) {
                    for (std::size_t j = 0; j < p; ++j) {
                        if (tp[h][j] < 0) continue;
                        std::vector<int> mm = m;
                        ++mm[tp[h][j]];
                        UP[k](a, pos[k + 1][mm] * p + j) += V(h, j);
                    }
                }
            }
        }
    }
    // Full bank: a completion hands the server to the next job, which starts in h.
    Matrix<T> CDEP(nk[cmax], nk[cmax], zero);
    for (std::size_t row = 0; row < nCfg[cmax]; ++row) {
        const std::vector<int>& m = cfg[cmax][row];
        for (std::size_t h = 0; h < p; ++h)
            for (std::size_t t = 0; t < NT; ++t) {
                if (m[t] == 0) continue;
                const T mult = num_traits<T>::from_int(m[t]);
                for (std::size_t j = 0; j < p; ++j) {
                    if (tp[h][j] < 0) continue;
                    std::vector<int> mm = m;
                    --mm[t];
                    ++mm[tp[h][j]];
                    CDEP(row * p + h, pos[cmax][mm] * p + j) += T(mult * done[t] * V(h, j));
                }
            }
    }
    std::vector<T> speed(Nlev + 1, one);
    if (!sf.empty()) {
        for (std::size_t n = 1; n <= Nlev; ++n) {
            const std::size_t b = std::min(n, cmax);
            const T bt = num_traits<T>::from_int(static_cast<int>(b));
            if (!(sf[n - 1] == bt)) speed[n] = T(sf[n - 1] / bt);
        }
    }
    LdqbdMphcBlocks<T> out;
    out.Q0.assign(Nlev, Matrix<T>(1, 1, zero));
    out.Q1.assign(Nlev + 1, Matrix<T>(1, 1, zero));
    out.Q2.assign(Nlev + 1, Matrix<T>(1, 1, zero));
    out.Q1[0] = Matrix<T>(p, p, zero);  // level 0: arrivals only, h frozen
    for (std::size_t i = 0; i < p; ++i) out.Q1[0](i, i) = T(-arrRate[0]);
    for (std::size_t n = 1; n <= Nlev; ++n) {
        const std::size_t b = std::min(n, cmax);
        Matrix<T> B(nk[b], nk[b], zero);
        for (std::size_t i = 0; i < nk[b]; ++i) {
            for (std::size_t j = 0; j < nk[b]; ++j) B(i, j) = T(speed[n] * LOC[b](i, j));
            B(i, i) -= arrRate[n];
        }
        out.Q1[n] = B;
    }
    for (std::size_t n = 0; n + 1 <= Nlev; ++n) {
        if (n < cmax) {
            Matrix<T> B(nk[n], nk[n + 1], zero);
            for (std::size_t i = 0; i < nk[n]; ++i)
                for (std::size_t j = 0; j < nk[n + 1]; ++j) B(i, j) = T(arrRate[n] * UP[n](i, j));
            out.Q0[n] = B;
        } else {
            Matrix<T> B(nk[cmax], nk[cmax], zero);
            for (std::size_t i = 0; i < nk[cmax]; ++i) B(i, i) = arrRate[n];
            out.Q0[n] = B;
        }
    }
    for (std::size_t n = 1; n <= Nlev; ++n) {
        const Matrix<T>& S = (n <= cmax) ? DN[n] : CDEP;
        Matrix<T> B(S.rows(), S.cols(), zero);
        for (std::size_t i = 0; i < S.rows(); ++i)
            for (std::size_t j = 0; j < S.cols(); ++j) B(i, j) = T(speed[n] * S(i, j));
        out.Q2[n] = B;
    }
    return out;
}

/**
 * Block-tridiagonal generator of an M/PH/c queue with level-dependent arrivals.
 *
 * @param D0      service sub-generator (p x p), phase changes without completion
 * @param D1      service completion block (p x p); D1 = (-D0*1)*alpha for a PH
 * @param alpha   length-p vector a server starts each new job in
 * @param c       number of identical servers (>= 1; capped at the top level)
 * @param arrRate length Nlev+1; arrRate[n] is the arrival rate out of level n
 * @param sf      empty, or length >= Nlev: a multiplier on the station's TOTAL
 *                service rate at level n (load dependence), so each busy server
 *                runs at sf[n-1]/min(n,c) of nominal and sf[n-1] = min(n,c)
 *                reproduces the unscaled queue exactly
 * @param carry   true for a correlated (non-renewal) MAP service, carried across services in
 *                start order, idle periods included (how JMT and LDES sample it); alpha is
 *                then unused. See `ldqbd_mapc_carry`
 *
 * Q2 carries an unused entry at index 0 so the three lists line up by level,
 * matching what the C++ `ldqbd` expects.
 */
template <class T>
LdqbdMphcBlocks<T> ldqbd_mphc(const Matrix<T>& D0, const Matrix<T>& D1,
                              const std::vector<T>& alpha, double c,
                              const std::vector<T>& arrRate, const std::vector<T>& sf,
                              bool carry = false) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t p = D0.rows();
    if (arrRate.empty())
        throw InputError("ldqbd_mphc: needs at least one level above the empty one");
    const std::size_t Nlev = arrRate.size() - 1;
    if (Nlev < 1)
        throw InputError("ldqbd_mphc: needs at least one level above the empty one");
    if (D0.cols() != p || D1.rows() != p || D1.cols() != p || alpha.size() != p)
        throw InputError("ldqbd_mphc: D0, D1 and alpha must all have the same order");
    if (!sf.empty() && sf.size() < Nlev)
        throw InputError("ldqbd_mphc: sf must give one total-service-rate factor per level");

    // Servers that can never be busy do not need a coordinate: above the top
    // level there is nothing left to serve.
    const std::size_t cmax =
        !std::isfinite(c) ? Nlev
                          : std::min(std::max<std::size_t>(1, static_cast<std::size_t>(
                                                                  std::llround(c))),
                                     Nlev);
    if (carry) return ldqbd_mapc_carry(D0, D1, cmax, arrRate, sf);

    std::vector<std::vector<std::vector<int> > > cfg(cmax + 1);
    std::vector<std::map<std::vector<int>, std::size_t> > pos(cmax + 1);
    std::vector<std::size_t> nCfg(cmax + 1, 0);
    for (std::size_t k = 0; k <= cmax; ++k) {
        cfg[k] = ph_multisets(p, k);
        for (std::size_t r = 0; r < cfg[k].size(); ++r) pos[k][cfg[k][r]] = r;
        nCfg[k] = cfg[k].size();
    }

    if (nCfg[cmax] > LDQBD_MPHC_MAX_CONFIGS)
        throw UnsupportedError(
            "ldqbd_mphc: the exact M/PH/c chain needs " + std::to_string(nCfg[cmax]) +
            " configurations per level for " + std::to_string(cmax) + " servers and " +
            std::to_string(p) +
            " service phases, above the " + std::to_string(LDQBD_MPHC_MAX_CONFIGS) +
            " the level-by-level inverses can carry. Use fewer phases (a lower-order fit), "
            "fewer servers, or SolverCTMC/SolverLDES on this model");

    // Completion rate out of each phase, summed over targets.
    std::vector<T> t(p, zero);
    for (std::size_t i = 0; i < p; ++i)
        for (std::size_t j = 0; j < p; ++j) t[i] += D1(i, j);

    // Structural blocks per busy-server count: LOC the within-level phase
    // changes (with the full outflow on its diagonal), UP the entry of a newly
    // busy server, DN a completion that leaves a server idle.
    std::vector<Matrix<T> > LOC(cmax + 1), UP(cmax + 1), DN(cmax + 1);
    for (std::size_t k = 0; k <= cmax; ++k) {
        const std::vector<std::vector<int> >& Ck = cfg[k];
        LOC[k] = Matrix<T>(nCfg[k], nCfg[k], zero);
        for (std::size_t row = 0; row < nCfg[k]; ++row) {
            const std::vector<int>& m = Ck[row];
            for (std::size_t i = 0; i < p; ++i) {
                if (m[i] == 0) continue;
                const T mult = num_traits<T>::from_int(m[i]);
                for (std::size_t j = 0; j < p; ++j) {
                    if (j == i) continue;
                    std::vector<int> mm = m;
                    --mm[i];
                    ++mm[j];
                    LOC[k](row, pos[k][mm]) += T(mult * D0(i, j));
                }
                // D0(i,i) is the total outflow of phase i, completions included
                LOC[k](row, row) += T(mult * D0(i, i));
            }
        }

        if (k < cmax) {
            UP[k] = Matrix<T>(nCfg[k], nCfg[k + 1], zero);
            for (std::size_t row = 0; row < nCfg[k]; ++row) {
                const std::vector<int>& m = Ck[row];
                for (std::size_t j = 0; j < p; ++j) {
                    std::vector<int> mm = m;
                    ++mm[j];
                    UP[k](row, pos[k + 1][mm]) += alpha[j];
                }
            }
        }

        if (k > 0) {
            DN[k] = Matrix<T>(nCfg[k], nCfg[k - 1], zero);
            for (std::size_t row = 0; row < nCfg[k]; ++row) {
                const std::vector<int>& m = Ck[row];
                for (std::size_t i = 0; i < p; ++i) {
                    if (m[i] == 0) continue;
                    std::vector<int> mm = m;
                    --mm[i];
                    DN[k](row, pos[k - 1][mm]) += T(num_traits<T>::from_int(m[i]) * t[i]);
                }
            }
        }
    }

    // A completion at a full server bank takes the next waiting job at once, so
    // the server stays busy and only its phase moves: the repeating down block.
    Matrix<T> CDEP(nCfg[cmax], nCfg[cmax], zero);
    for (std::size_t row = 0; row < nCfg[cmax]; ++row) {
        const std::vector<int>& m = cfg[cmax][row];
        for (std::size_t i = 0; i < p; ++i) {
            if (m[i] == 0) continue;
            const T mult = num_traits<T>::from_int(m[i]);
            for (std::size_t j = 0; j < p; ++j) {
                std::vector<int> mm = m;
                --mm[i];
                ++mm[j];
                CDEP(row, pos[cmax][mm]) += T(mult * D1(i, j));
            }
        }
    }

    // Per-server speed. Without load dependence every busy server runs at its
    // nominal rate; with it, the aggregate sf(n) is shared over the busy
    // servers. sf(n) == min(n,c) is passed through as exactly one so the
    // unscaled chain is reproduced bit for bit.
    const T one = num_traits<T>::from_int(1);
    std::vector<T> speed(Nlev + 1, one);
    if (!sf.empty()) {
        for (std::size_t n = 1; n <= Nlev; ++n) {
            const std::size_t b = std::min(n, cmax);
            const T bt = num_traits<T>::from_int(static_cast<int>(b));
            if (!(sf[n - 1] == bt)) speed[n] = T(sf[n - 1] / bt);
        }
    }

    LdqbdMphcBlocks<T> out;
    out.Q0.assign(Nlev, Matrix<T>(1, 1, zero));
    out.Q1.assign(Nlev + 1, Matrix<T>(1, 1, zero));
    out.Q2.assign(Nlev + 1, Matrix<T>(1, 1, zero));

    out.Q1[0] = Matrix<T>(1, 1, T(-arrRate[0]));  // level 0: arrivals only
    for (std::size_t n = 1; n <= Nlev; ++n) {
        const std::size_t b = std::min(n, cmax);
        Matrix<T> B(nCfg[b], nCfg[b], zero);
        for (std::size_t i = 0; i < nCfg[b]; ++i) {
            for (std::size_t j = 0; j < nCfg[b]; ++j) B(i, j) = T(speed[n] * LOC[b](i, j));
            B(i, i) -= arrRate[n];
        }
        out.Q1[n] = B;
    }

    for (std::size_t n = 0; n + 1 <= Nlev; ++n) {
        if (n < cmax) {
            // a free server takes the job, starting it in a phase drawn from alpha
            Matrix<T> B(nCfg[n], nCfg[n + 1], zero);
            for (std::size_t i = 0; i < nCfg[n]; ++i)
                for (std::size_t j = 0; j < nCfg[n + 1]; ++j) B(i, j) = T(arrRate[n] * UP[n](i, j));
            out.Q0[n] = B;
        } else {
            // the job waits, so every busy phase is unchanged
            Matrix<T> B(nCfg[cmax], nCfg[cmax], zero);
            for (std::size_t i = 0; i < nCfg[cmax]; ++i) B(i, i) = arrRate[n];
            out.Q0[n] = B;
        }
    }

    for (std::size_t n = 1; n <= Nlev; ++n) {
        if (n <= cmax) {
            Matrix<T> B(nCfg[n], nCfg[n - 1], zero);   // the server falls idle
            for (std::size_t i = 0; i < nCfg[n]; ++i)
                for (std::size_t j = 0; j < nCfg[n - 1]; ++j) B(i, j) = T(speed[n] * DN[n](i, j));
            out.Q2[n] = B;
        } else {
            Matrix<T> B(nCfg[cmax], nCfg[cmax], zero);  // it takes the next job
            for (std::size_t i = 0; i < nCfg[cmax]; ++i)
                for (std::size_t j = 0; j < nCfg[cmax]; ++j) B(i, j) = T(speed[n] * CDEP(i, j));
            out.Q2[n] = B;
        }
    }

    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_LDQBD_MPHC_H
