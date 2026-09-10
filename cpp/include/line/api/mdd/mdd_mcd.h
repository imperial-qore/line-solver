/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_MCD_H
#define LINE_API_MDD_MDD_MCD_H

/**
 * Miner-Ciardo-Donatelli approximate stationary analysis.
 *
 * Port of matlab/src/api/mdd/mdd_mcd.m, jline.api.mdd.Mdd_mcd and
 * python/line_solver/api/mdd/mcd.py, after A.S. Miner, G. Ciardo, S. Donatelli,
 * "Using the exact state space of a Markov model to compute approximate
 * stationary measures", ACM SIGMETRICS 2000, pp.207-216.
 *
 * Solve a structured CTMC whose EXACT reachable state space is stored in a
 * decision diagram, by building and iterating K level-CTMCs. The method never
 * forms the |S|-state generator or probability vector. It keeps one CTMC per
 * level k, over states M_k = {(p,i_k)} with p a level-k node and i_k a local
 * state on a non-null arc, and iterates the coupled system to a fixed point.
 * The single approximation (Eq. 5) is Pr{i_k | alpha} = Pr{i_k | p}: the
 * local-state law at level k depends only on the node p, not the full path
 * above it, which the exact reachability the diagram encodes justifies. For
 * product-form models the method is EXACT (paper Sec. 5), so on a single-class
 * closed QN it reproduces SolverCTMC.
 *
 * ORIENTATION. The paper indexes levels K (top/root) down to 1
 * (bottom/terminal); `MDD` uses level 0 as the root. This function works in the
 * paper's orientation with 0-based indices, so paper level k (0 = bottom) maps
 * to MDD level K-1-k and to station K-1-k. Getting this backwards silently
 * mislabels every per-station metric.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lstsq.h"
#include "line/util/matrix.h"

namespace line {
namespace mdd {

namespace detail {

/** Node marginal of a level vector: Pr{p} = sum over the arcs of p. */
template <class T>
std::vector<T> mcd_node_marginal(const std::vector<std::pair<int, int>>& rows,
                                 const std::vector<T>& pk, int nnodes) {
    std::vector<T> pr(static_cast<std::size_t>(nnodes), num_traits<T>::from_int(0));
    for (std::size_t r = 0; r < rows.size(); ++r) pr[rows[r].first - 1] += pk[r];
    return pr;
}

template <class T>
std::vector<std::vector<T>> mcd_identity(std::size_t n) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<std::vector<T>> I(n, std::vector<T>(n, zero));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = one;
    return I;
}

/** Infinitesimal generator of a rate matrix: the diagonal absorbs the row sum. */
template <class T>
std::vector<std::vector<T>> mcd_generator(const std::vector<std::vector<T>>& R) {
    const std::size_t n = R.size();
    const T zero = num_traits<T>::from_int(0);
    std::vector<std::vector<T>> Q = R;
    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += R[i][j];
        Q[i][i] -= s;
    }
    return Q;
}

/**
 * Least-squares solution of an overdetermined system by Householder QR.
 *
 * QR is used rather than the normal equations because A'A squares the condition
 * number, which is exactly the failure the appended-normalisation form of
 * `mcd_solve_stat` exists to avoid. This is why the shared `line::lstsq` is not
 * called here: it forms the normal equations on the full-rank branch.
 */
template <class T>
std::vector<T> mcd_lstsq(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
    using std::sqrt;
    const std::size_t m = A.size(), n = A[0].size();
    const T zero = num_traits<T>::from_int(0), two = num_traits<T>::from_int(2);
    std::vector<std::vector<T>> R = A;
    std::vector<T> y = b;
    for (std::size_t k = 0; k < n; ++k) {
        T norm2 = zero;
        for (std::size_t i = k; i < m; ++i) norm2 += T(R[i][k] * R[i][k]);
        if (norm2 == zero) continue;
        T norm = T(sqrt(norm2));
        if (R[k][k] > zero) norm = T(-norm);
        std::vector<T> v(m, zero);
        for (std::size_t i = k; i < m; ++i) v[i] = R[i][k];
        v[k] -= norm;
        T vtv = zero;
        for (std::size_t i = k; i < m; ++i) vtv += T(v[i] * v[i]);
        if (vtv == zero) continue;
        for (std::size_t j = k; j < n; ++j) {
            T dot = zero;
            for (std::size_t i = k; i < m; ++i) dot += T(v[i] * R[i][j]);
            const T f = T(two * dot / vtv);
            for (std::size_t i = k; i < m; ++i) R[i][j] -= T(f * v[i]);
        }
        T dot = zero;
        for (std::size_t i = k; i < m; ++i) dot += T(v[i] * y[i]);
        const T f = T(two * dot / vtv);
        for (std::size_t i = k; i < m; ++i) y[i] -= T(f * v[i]);
    }
    std::vector<T> x(n, zero);
    for (std::size_t ii = n; ii > 0; --ii) {
        const std::size_t i = ii - 1;
        T s = y[i];
        for (std::size_t j = i + 1; j < n; ++j) s -= T(R[i][j] * x[j]);
        x[i] = R[i][i] == zero ? zero : T(s / R[i][i]);
    }
    return x;
}

/**
 * Stationary distribution of a small irreducible generator: p Q = 0, sum p = 1.
 *
 * The normalisation is APPENDED rather than substituted for the last balance
 * equation: overwriting a row discards a constraint and leaves Q' singular to
 * working precision from about |M_k| = 325 upwards, so the solve returns NaN.
 * The overdetermined system has full column rank whenever the level chain is
 * irreducible, and least squares solves it stably.
 *
 * TWO BACKENDS, chosen by the arithmetic rather than by an option. In floating
 * point the Householder QR above is used, because the shared `line::lstsq`
 * forms the normal equations on its full-rank branch and A'A squares the
 * condition number -- exactly the failure the appended normalisation exists to
 * avoid. Under EXACT arithmetic there is no condition number to square and no
 * square root to take, so `line::lstsq` is called directly: that is what makes
 * the level aggregation available at `Rational`, where it returns the fixed
 * point of the level system with no rounding at all.
 */
template <class T>
std::vector<T> mcd_solve_stat(const std::vector<std::vector<T>>& Q) {
    const std::size_t n = Q.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (n == 1) return std::vector<T>(1, one);
    std::vector<std::vector<T>> A(n + 1, std::vector<T>(n, zero));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A[j][i] = Q[i][j];
    for (std::size_t j = 0; j < n; ++j) A[n][j] = one;
    std::vector<T> rhs(n + 1, zero);
    rhs[n] = one;
    std::vector<T> p;
    if constexpr (num_traits<T>::is_exact) {
        Matrix<T> Am(n + 1, n, zero);
        for (std::size_t i = 0; i <= n; ++i)
            for (std::size_t j = 0; j < n; ++j) Am(i, j) = A[i][j];
        p = line::lstsq(Am, rhs).x;
    } else {
        p = mcd_lstsq(A, rhs);
    }
    T s = zero;
    for (std::size_t i = 0; i < n; ++i) {
        if (p[i] < zero) p[i] = zero;
        s += p[i];
    }
    const double sd = num_traits<T>::to_double(s);
    if (!(sd > 0) || std::isnan(sd) || std::isinf(sd))
        throw NumericError("mdd_mcd: level CTMC of order " + std::to_string(n) +
                           " admits no proper stationary distribution (the level generator is "
                           "reducible or numerically degenerate).");
    for (std::size_t i = 0; i < n; ++i) p[i] = T(p[i] / s);
    return p;
}

/**
 * Path counts per MDD level: above[oL][p] is the number of distinct root-to-p
 * paths, |A(p)| in the paper's notation, and below[oL][p] the number of
 * accepted states under p.
 *
 * Both are O(#nodes) and serve two purposes: the uniform initialisation, and
 * the exactness certificate. The single approximation is
 * Pr{i_k | alpha} = Pr{i_k | p}, so when a node is reached by exactly one path,
 * conditioning on the node IS conditioning on the path and the identity is
 * exact. That test is SUFFICIENT, not necessary: a product-form model is exact
 * too however much its diagram shares. Note also that max |A(p)| = 1 means no
 * node is shared, i.e. the diagram compresses nothing, so exactness by this
 * route and a useful saving are mutually exclusive.
 */
inline void mcd_path_counts(const MddStruct& mdds, std::size_t K,
                            std::vector<std::vector<double>>& above,
                            std::vector<std::vector<double>>& below) {
    below.assign(K, std::vector<double>());
    above.assign(K, std::vector<double>());
    for (std::size_t oo = K; oo > 0; --oo) {
        const std::size_t oL = oo - 1;
        std::vector<double> nb(static_cast<std::size_t>(mdds.nnodes[oL]), 0.0);
        for (std::size_t p = 0; p < nb.size(); ++p) {
            double s = 0;
            for (int v = 0; v < mdds.domain[oL]; ++v) {
                const int ch = mdds.node[oL][p][v];
                if (oL + 1 == K) {
                    if (ch == TERM_TRUE) s += 1.0;
                } else if (ch > 0) {
                    s += below[oL + 1][ch - 1];
                }
            }
            nb[p] = s;
        }
        below[oL] = nb;
    }
    for (std::size_t oL = 0; oL < K; ++oL)
        above[oL].assign(static_cast<std::size_t>(mdds.nnodes[oL]), 0.0);
    above[0][mdds.root - 1] = 1.0;
    for (std::size_t oL = 0; oL + 1 < K; ++oL) {
        for (std::size_t p = 0; p < above[oL].size(); ++p) {
            const double w = above[oL][p];
            if (w == 0) continue;
            for (int v = 0; v < mdds.domain[oL]; ++v) {
                const int ch = mdds.node[oL][p][v];
                if (ch > 0) above[oL + 1][ch - 1] += w;
            }
        }
    }
}

/**
 * Uniform law over the EXACT reachable set, projected onto each level.
 *
 * Pr{(p,v)} = (paths root->p) * (states below arc p[v]) / |S|. A flat law over
 * M_k instead treats level states as equiprobable irrespective of how many
 * global states they stand for, which breaks the population invariant the
 * diagram encodes; from about K=8 the coupled iteration then descends into the
 * basin of the DEGENERATE empty-population fixed point (all mass on local state
 * 0 at every level, a true fixed point since no station can then emit) and
 * converges to it with zero residual. The projection below is consistent across
 * levels by construction, so the iteration starts inside the physical simplex.
 */
template <class T>
std::vector<std::vector<T>> mcd_uniform_init(
    const MddStruct& mdds, const std::vector<std::vector<std::pair<int, int>>>& Mrows,
    std::size_t K, const std::vector<std::vector<double>>& above,
    const std::vector<std::vector<double>>& below) {
    std::vector<std::vector<T>> pik(K);
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t oL = K - 1 - k;  // paper level k is MDD level K-1-k
        const std::vector<std::pair<int, int>>& rows = Mrows[k];
        std::vector<double> w(rows.size(), 0.0);
        double sum = 0;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const int p = rows[r].first, v = rows[r].second;
            double val;
            if (oL + 1 == K) {
                val = above[oL][p - 1];  // a TRUE arc stands for one state
            } else {
                const int ch = mdds.node[oL][p - 1][v];
                val = above[oL][p - 1] * below[oL + 1][ch - 1];
            }
            w[r] = val;
            sum += val;
        }
        pik[k].assign(rows.size(), num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < rows.size(); ++r)
            pik[k][r] = num_traits<T>::from_double(w[r] / sum);
    }
    return pik;
}

/**
 * ComputeAs(k): A_k^e from A_{k+1}^e, the "from above" contribution (Fig. 3).
 *
 * The adjust denominator Pr{p[v]} is the FROM-ABOVE marginal of the child node,
 * Pr{p} = sum over parents of pi_{k+1}. Using it (rather than the level-k CTMC
 * marginal, which only equals it at convergence) makes adjust a proper
 * conditional Pr{(parent,arc)|child} and pins the inter-level node marginals,
 * removing the spurious fixed points.
 */
template <class T>
std::vector<std::vector<std::vector<T>>> mcd_compute_as(
    std::size_t k, const std::vector<std::vector<std::vector<T>>>& Aup,
    const std::vector<std::vector<std::vector<int>>>& Pnode,
    const std::vector<std::vector<T>>& pik, const std::vector<std::vector<MddLocalMatrix<T>>>& W,
    const std::vector<std::vector<std::pair<int, int>>>& Mrows, const std::vector<int>& nn,
    std::size_t E) {
    const T zero = num_traits<T>::from_int(0);
    const std::vector<std::pair<int, int>>& rows1 = Mrows[k + 1];
    std::vector<T> pr_above(static_cast<std::size_t>(nn[k]), zero);
    for (std::size_t r = 0; r < rows1.size(); ++r) {
        const int child = Pnode[k + 1][rows1[r].first - 1][rows1[r].second];
        if (child > 0) pr_above[child - 1] += pik[k + 1][r];
    }
    std::vector<std::vector<std::vector<T>>> Ak(
        E, std::vector<std::vector<T>>(static_cast<std::size_t>(nn[k]),
                                       std::vector<T>(static_cast<std::size_t>(nn[k]), zero)));
    for (std::size_t r = 0; r < rows1.size(); ++r) {
        const int p = rows1[r].first;
        const int v = rows1[r].second;
        const int childp = Pnode[k + 1][p - 1][v];  // p[v]: node at level k
        if (childp <= 0 || !(pr_above[childp - 1] > zero)) continue;
        const T adjust = T(pik[k + 1][r] / pr_above[childp - 1]);
        for (std::size_t e = 0; e < E; ++e) {
            const std::vector<std::size_t>& wcols = W[e][k + 1].cols[v];
            if (wcols.empty()) continue;
            const std::vector<T>& wvals = W[e][k + 1].vals[v];
            const std::vector<T>& arow = Aup[e][p - 1];
            for (std::size_t wi = 0; wi < wcols.size(); ++wi) {
                const std::size_t w = wcols[wi];
                const T wv = wvals[wi];
                for (std::size_t q = 0; q < arow.size(); ++q) {
                    if (arow[q] == zero) continue;
                    const int childq = Pnode[k + 1][q][w];
                    if (childq <= 0) continue;  // q[w] null
                    Ak[e][childp - 1][childq - 1] += T(arow[q] * wv * adjust);
                }
            }
        }
    }
    return Ak;
}

/**
 * ComputeMC(k): level-k rate matrix (Eq. 6),
 * R_k^e[(p,i),(q,j)] = A_k^e[p,q] * W_k^e[i,j] * b_{k-1}^e[p[i]].
 */
template <class T>
std::vector<std::vector<T>> mcd_compute_mc(
    std::size_t k, const std::vector<std::vector<std::vector<T>>>& Ak,
    const std::vector<std::vector<std::vector<T>>>& bcell,
    const std::vector<std::vector<std::vector<int>>>& Pnode,
    const std::vector<std::vector<MddLocalMatrix<T>>>& W,
    const std::vector<std::vector<std::pair<int, int>>>& Mrows,
    const std::vector<std::vector<int>>& Midx, const std::vector<std::size_t>& level_sizes,
    const std::vector<int>& dom, std::size_t E) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t nm = level_sizes[k];
    const std::vector<std::pair<int, int>>& rows = Mrows[k];
    std::vector<std::vector<T>> Rk(nm, std::vector<T>(nm, zero));
    for (std::size_t r = 0; r < nm; ++r) {
        const int p = rows[r].first;
        const int v = rows[r].second;
        for (std::size_t e = 0; e < E; ++e) {
            const std::vector<std::size_t>& wcols = W[e][k].cols[v];
            if (wcols.empty()) continue;
            T bfac = one;  // terminal ONE
            if (k > 0) bfac = bcell[k - 1][Pnode[k][p - 1][v] - 1][e];
            if (bfac == zero) continue;
            const std::vector<T>& wvals = W[e][k].vals[v];
            const std::vector<T>& arow = Ak[e][p - 1];
            for (std::size_t wi = 0; wi < wcols.size(); ++wi) {
                const std::size_t w = wcols[wi];
                const T wv = wvals[wi];
                for (std::size_t q = 0; q < arow.size(); ++q) {
                    if (arow[q] == zero) continue;
                    const int di = Midx[k][q * static_cast<std::size_t>(dom[k]) + w];
                    if (di == 0) continue;  // (q,w) not in M_k
                    Rk[r][di - 1] += T(arow[q] * wv * bfac);
                }
            }
        }
    }
    return Rk;
}

}  // namespace detail

/**
 * Approximate stationary measures by decision-diagram-guided aggregation.
 *
 * @param mdds the reachable set, in MDD orientation
 * @param desc the Kronecker rate descriptor
 * @param options level-iteration knobs
 */
template <class T>
MddMcdResult<T> mdd_mcd(const MddStruct& mdds, const MddDescriptor<T>& desc,
                        const MddMcdOptions& options = MddMcdOptions()) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = mdds.K;
    if (K == 0) throw InputError("mdd_mcd: the diagram has no levels");

    // ---- paper orientation: paper level k <-> MDD level K-1-k = station K-1-k
    std::vector<std::vector<std::vector<int>>> Pnode(K);
    std::vector<int> nn(K, 0), dom(K, 0);
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t oL = K - 1 - k;
        Pnode[k] = mdds.node[oL];
        nn[k] = mdds.nnodes[oL];
        dom[k] = mdds.domain[oL];
    }

    // ---- per (event, paper level) local matrices W_k^e; an untouched level
    // carries the identity, which is supplied here rather than stored
    const std::size_t E = desc.events.size();
    std::vector<std::vector<MddLocalMatrix<T>>> W(E, std::vector<MddLocalMatrix<T>>(K));
    std::vector<std::vector<bool>> touched(E, std::vector<bool>(K, false));
    for (std::size_t e = 0; e < E; ++e) {
        const MddEvent<T>& ev = desc.events[e];
        for (std::size_t t = 0; t < ev.lev.size(); ++t) {
            const std::size_t pk = K - 1 - ev.lev[t];  // station level -> paper level
            W[e][pk] = ev.W[t];
            touched[e][pk] = true;
        }
        for (std::size_t k = 0; k < K; ++k)
            if (!touched[e][k])
                W[e][k] = MddLocalMatrix<T>::identity(static_cast<std::size_t>(dom[k]));
    }

    // ---- level-k CTMC state sets M_k = {(p, v) : arc p[v] non-null}
    std::vector<std::vector<std::pair<int, int>>> Mrows(K);
    std::vector<std::vector<int>> Midx(K);
    std::vector<std::size_t> level_sizes(K, 0);
    for (std::size_t k = 0; k < K; ++k) {
        // column-major order, matching the MATLAB find() the reference uses
        for (int v = 0; v < dom[k]; ++v) {
            for (int p = 0; p < nn[k]; ++p) {
                const int ch = Pnode[k][p][v];
                const bool live = (k == 0) ? (ch == TERM_TRUE) : (ch > 0);
                if (live) Mrows[k].push_back(std::make_pair(p + 1, v));
            }
        }
        level_sizes[k] = Mrows[k].size();
        Midx[k].assign(static_cast<std::size_t>(nn[k]) * static_cast<std::size_t>(dom[k]), 0);
        for (std::size_t r = 0; r < Mrows[k].size(); ++r)
            Midx[k][static_cast<std::size_t>(Mrows[k][r].first - 1) *
                        static_cast<std::size_t>(dom[k]) +
                    static_cast<std::size_t>(Mrows[k][r].second)] = static_cast<int>(r + 1);
    }

    // ---- initialise level stationary vectors and node marginals
    std::vector<std::vector<double>> above, below;
    detail::mcd_path_counts(mdds, K, above, below);
    std::vector<std::vector<T>> pik;
    if (!options.initpik.empty()) {
        pik.assign(K, std::vector<T>());
        for (std::size_t k = 0; k < K; ++k) {
            pik[k].assign(options.initpik[k].size(), zero);
            for (std::size_t r = 0; r < options.initpik[k].size(); ++r)
                pik[k][r] = num_traits<T>::from_double(options.initpik[k][r]);
        }
    } else {
        pik = detail::mcd_uniform_init<T>(mdds, Mrows, K, above, below);
    }
    std::vector<std::vector<T>> Prp(K);
    for (std::size_t k = 0; k < K; ++k) Prp[k] = detail::mcd_node_marginal(Mrows[k], pik[k], nn[k]);

    // ---- fixed-point iteration (Fig. 3, procedure Solve)
    int iters = 0;
    bool converged = false;
    double delta = 0;
    for (int it = 1; it <= options.maxiter; ++it) {
        iters = it;
        const std::vector<std::vector<T>> piold = pik;

        // ComputeBs, bottom-up:
        // b_k^e[p] = sum_v Pr{v|p} * b_{k-1}^e[p[v]] * lambda
        std::vector<std::vector<std::vector<T>>> bcell(K);
        for (std::size_t k = 0; k < K; ++k) {
            std::vector<std::vector<T>> bk(static_cast<std::size_t>(nn[k]),
                                           std::vector<T>(E, zero));
            const std::vector<std::pair<int, int>>& rows = Mrows[k];
            for (std::size_t r = 0; r < rows.size(); ++r) {
                const int p = rows[r].first;
                const int v = rows[r].second;
                if (!(Prp[k][p - 1] > zero)) continue;
                const T adjust = T(pik[k][r] / Prp[k][p - 1]);  // Pr{v|p}
                for (std::size_t e = 0; e < E; ++e) {
                    const T le = W[e][k].row_sum[v];
                    if (le == zero) continue;  // not locally enabled
                    T down = one;              // terminal ONE
                    if (k > 0) down = bcell[k - 1][Pnode[k][p - 1][v] - 1][e];
                    bk[p - 1][e] += T(adjust * down * le);
                }
            }
            bcell[k] = bk;
        }

        // top-down: ComputeAs(k) then SolveLevel(k)
        std::vector<std::vector<std::vector<std::vector<T>>>> Acell(K);
        Acell[K - 1].assign(E, std::vector<std::vector<T>>());
        for (std::size_t e = 0; e < E; ++e)
            Acell[K - 1][e] = detail::mcd_identity<T>(static_cast<std::size_t>(nn[K - 1]));
        for (std::size_t kk = K; kk > 0; --kk) {
            const std::size_t k = kk - 1;
            if (k + 1 < K)
                Acell[k] = detail::mcd_compute_as(k, Acell[k + 1], Pnode, pik, W, Mrows, nn, E);
            // SolveLevel(k): assemble R_k (Eq. 6), solve pi_k Q_k = 0
            const std::vector<std::vector<T>> Rk =
                detail::mcd_compute_mc(k, Acell[k], bcell, Pnode, W, Mrows, Midx, level_sizes,
                                       dom, E);
            pik[k] = detail::mcd_solve_stat(detail::mcd_generator(Rk));
            Prp[k] = detail::mcd_node_marginal(Mrows[k], pik[k], nn[k]);
        }

        // A diverged iterate must not be read as converged, which is what a
        // NaN-skipping maximum would do.
        delta = 0;
        for (std::size_t k = 0; k < K; ++k) {
            double dk = 0;
            for (std::size_t r = 0; r < pik[k].size(); ++r) {
                const double d = num_traits<T>::to_double(T(pik[k][r] - piold[k][r]));
                dk = std::max(dk, d < 0 ? -d : d);
            }
            if (std::isnan(dk) || std::isinf(dk))
                throw NumericError("mdd_mcd: level " + std::to_string(k + 1) +
                                   " iterate is not finite at iteration " + std::to_string(it) +
                                   "; the level CTMC did not yield a proper stationary vector.");
            delta = std::max(delta, dk);
        }
        if (delta < options.tol) {
            converged = true;
            break;
        }
    }
    if (!converged)
        throw NumericError("mdd_mcd: the coupled level iteration did not converge in " +
                           std::to_string(options.maxiter) + " sweeps (last change " +
                           std::to_string(delta) + " against tol " + std::to_string(options.tol) +
                           "); the level marginals returned would not be a fixed point. Raise "
                           "maxiter or relax tol.");

    // ---- performance measures from the per-level marginals. QLen is the mean
    // local value and is defined for any descriptor (jobs at a station, tokens
    // in a place); X and U need the queueing parameters.
    const bool is_qn = !desc.mu.empty() && !desc.servers.empty();
    MddMcdResult<T> out;
    out.QLen.assign(K, zero);
    if (is_qn) {
        out.X.assign(K, zero);
        out.U.assign(K, zero);
    }
    for (std::size_t s = 0; s < K; ++s) {
        const std::size_t k = K - 1 - s;  // paper level of station s
        const std::vector<std::pair<int, int>>& rows = Mrows[k];
        T q = zero, busy = zero;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const double vd = desc.valuemap.empty()
                                  ? static_cast<double>(rows[r].second)
                                  : desc.valuemap[s][static_cast<std::size_t>(rows[r].second)];
            const T v = num_traits<T>::from_double(vd);
            q += T(v * pik[k][r]);
            if (is_qn) {
                const double srv = desc.servers[s];
                const T cap = num_traits<T>::from_double(vd < srv ? vd : srv);
                busy += T(cap * pik[k][r]);
            }
        }
        out.QLen[s] = q;
        if (is_qn) {
            out.X[s] = T(desc.mu[s] * busy);
            out.U[s] = std::isinf(desc.servers[s])
                           ? q
                           : T(busy / num_traits<T>::from_double(desc.servers[s]));
        }
    }

    // The level chains are coupled only through rates, so nothing in the
    // iteration forces the marginals to describe the same population; a fixed
    // point that does not is a wrong answer, not an approximation, and must not
    // be returned. The test is a conservation law of the model: the closed
    // population for a QN, a place invariant w'*m = const for a net.
    std::vector<double> winv = desc.invariant_weights;
    double vinv = desc.invariant_value;
    if (winv.empty() && desc.N > 0) {
        winv.assign(K, 1.0);  // closed QN: total population
        vinv = desc.N;
    }
    if (!winv.empty()) {
        double got = 0;
        for (std::size_t s = 0; s < K; ++s)
            got += winv[s] * num_traits<T>::to_double(out.QLen[s]);
        const double scale = std::max(1.0, std::fabs(vinv));
        if (std::fabs(got - vinv) > 1e-6 * scale)
            throw NumericError("mdd_mcd: the level marginals converged to an invariant value of " +
                               std::to_string(got) + " against the model value " +
                               std::to_string(vinv) + ", so the fixed point reached is degenerate "
                               "(the level chains are mutually inconsistent). Supply initpik with "
                               "a consistent starting law.");
    }

    out.pik = pik;
    out.Mrows = Mrows;
    out.level_sizes = level_sizes;
    out.iters = iters;
    out.paths_per_level.assign(K, 1.0);
    bool no_agg = true;
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t oL = K - 1 - k;
        double mx = 1.0;
        for (std::size_t p = 0; p < above[oL].size(); ++p) mx = std::max(mx, above[oL][p]);
        out.paths_per_level[k] = mx;
        if (mx > 1.0 + 1e-12) no_agg = false;
    }
    out.no_aggregation = no_agg;
    return out;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_MCD_H
