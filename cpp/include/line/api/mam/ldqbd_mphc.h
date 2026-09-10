/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_LDQBD_MPHC_H
#define LINE_API_MAM_LDQBD_MPHC_H

/**
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
#include <vector>

#include "line/util/error.h"
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
 *
 * Q2 carries an unused entry at index 0 so the three lists line up by level,
 * matching what the C++ `ldqbd` expects.
 */
template <class T>
LdqbdMphcBlocks<T> ldqbd_mphc(const Matrix<T>& D0, const Matrix<T>& D1,
                              const std::vector<T>& alpha, double c,
                              const std::vector<T>& arrRate, const std::vector<T>& sf) {
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
