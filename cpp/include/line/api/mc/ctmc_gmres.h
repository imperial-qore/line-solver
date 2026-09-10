/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_GMRES_H
#define LINE_API_MC_CTMC_GMRES_H

/**
 * Restarted GMRES with an ILUT preconditioner, for the linear systems a
 * generator produces.
 *
 * Templated port of matlab/src/api/mc/ctmc_gmres.m and
 * jar/src/main/java/jline/api/mc/Ctmc_gmres.java. This is the iterative
 * counterpart of the direct solve in ctmc_solve, meant for generators whose LU
 * fill-in exceeds the memory available; no CTMC-specific processing happens
 * here, so the same routine serves the stochastic complement and the
 * aggregation methods.
 *
 * Two preparation steps are not optional on a generator, and both are ported.
 * Rows are equilibrated to unit max norm, so the O(1) normalization row does
 * not mix with rows carrying rates of a wholly different magnitude. The states
 * are then reordered by reverse Cuthill-McKee: in the natural ordering of a
 * birth-death chain the unpivoted incomplete elimination has growth factor
 * (mu/lambda)^n, which overflows within a few thousand states, and a
 * fill-reducing ordering rather than pivoting is what removes it.
 *
 * The preconditioner is ILUT(p, tau) of Saad (1994): the row is expanded into a
 * dense workspace, entries below tau times the mean magnitude of the original
 * row are dropped as they are produced, and the p largest survivors are kept in
 * each of the L and U parts. A vanishing pivot is replaced by a value of the
 * order of the row threshold, keeping its sign, which is Saad's remedy and
 * avoids abandoning the factorization for one rate-free state. When the
 * factorization breaks down entirely the preconditioner degrades to Jacobi.
 *
 * MATLAB-VS-JAVA DISAGREEMENT, resolved in favour of Java. MATLAB calls the
 * built-in gmres with (L,U), which applies the preconditioner on the LEFT, so
 * its RELRES is the preconditioned residual norm(M\\(b-Ax))/norm(M\\b) and
 * depends on the preconditioner. The JAR expands the Krylov space of A M^-1
 * instead, and reports the true residual norm(b-Ax)/norm(b). The solution both
 * converge to is the same; the reported relres is not, and a tolerance on the
 * true residual is the one a caller can act on, so this port follows the JAR.
 * FLAG keeps the MATLAB convention: 0 converged, 1 iteration limit, 3
 * stagnation or divergence or a non-finite iterate.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC. GMRES stops on a residual tolerance and
 * its Arnoldi step normalizes by a Euclidean norm, so square roots appear in
 * every iteration and there is no exact answer to converge to: at Rational the
 * iteration would run to the iteration limit with exploding denominators. The
 * exact solve of the same system is ctmc_solve, which is what a Rational
 * caller should use.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <queue>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * Order above which the direct sparse factorization is abandoned in favour of
 * the Krylov path. The same number the other three codebases use
 * (`ctmc_solve.m`, `Ctmc_solve.GMRES_MIN_STATES`, `ctmc.py`), so a model
 * switches methods at the same size in all four.
 */
constexpr std::size_t GMRES_MIN_STATES = 6000;

template <class T>
struct GmresResult {
    std::vector<T> x;  ///< solution
    int flag;          ///< 0 converged, 1 iteration limit, 3 stagnation/divergence
    T relres;          ///< true relative residual norm(b - A x) / norm(b)
    long iter;         ///< total inner iterations
};

namespace detail {

/** Relative threshold below which ILUT discards a fill-in entry. */
constexpr double GMRES_ILUT_DROP_TOL = 1e-4;
/** Fill allowed per row, as a multiple of the average nnz per row. */
constexpr double GMRES_ILUT_FILL_FACTOR = 10.0;
/** Arnoldi breakdown threshold, relative to the unorthogonalized vector. */
constexpr double GMRES_BREAKDOWN_TOL = 1e-14;
constexpr double GMRES_DEFAULT_TOL = 1e-12;
constexpr long GMRES_DEFAULT_RESTART = 50;

/**
 * Finiteness without depending on which of std:: or boost:: supplies isfinite
 * for T: |v| < |v| + 1 holds for every finite value, fails for an infinity
 * (inf + 1 is inf) and fails for a NaN (every comparison does).
 */
template <class T>
inline bool num_isfinite(const T& v) {
    const T a = num_abs(T(v));
    return a < T(a + num_traits<T>::from_int(1));
}

/**
 * Compressed sparse row form, with each row's column indices increasing and the
 * position of the diagonal recorded. The incomplete factorization is a
 * row-oriented elimination and needs exactly this layout.
 */
template <class T>
struct CsrMatrix {
    std::size_t n = 0;
    std::vector<std::size_t> rowPtr, colIdx;
    std::vector<T> val;
    std::vector<long> diagPtr;

    void build_diag() {
        diagPtr.assign(n, -1);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t p = rowPtr[i]; p < rowPtr[i + 1]; ++p)
                if (colIdx[p] == i) {
                    diagPtr[i] = static_cast<long>(p);
                    break;
                }
    }

    static CsrMatrix of(const Matrix<T>& A) {
        CsrMatrix c;
        c.n = A.rows();
        const T zero = num_traits<T>::from_int(0);
        c.rowPtr.assign(c.n + 1, 0);
        for (std::size_t i = 0; i < c.n; ++i) {
            std::size_t cnt = 0;
            for (std::size_t j = 0; j < A.cols(); ++j)
                if (A(i, j) != zero) ++cnt;
            c.rowPtr[i + 1] = c.rowPtr[i] + cnt;
        }
        c.colIdx.resize(c.rowPtr[c.n]);
        c.val.resize(c.rowPtr[c.n]);
        std::size_t p = 0;
        for (std::size_t i = 0; i < c.n; ++i)
            for (std::size_t j = 0; j < A.cols(); ++j)
                if (A(i, j) != zero) {
                    c.colIdx[p] = j;
                    c.val[p] = A(i, j);
                    ++p;
                }
        c.build_diag();
        return c;
    }

    void mult(const std::vector<T>& v, std::vector<T>& out) const {
        for (std::size_t i = 0; i < n; ++i) {
            T s = num_traits<T>::from_int(0);
            for (std::size_t p = rowPtr[i]; p < rowPtr[i + 1]; ++p) s += val[p] * v[colIdx[p]];
            out[i] = s;
        }
    }

    /** Scales every row to unit max norm; returns the divisors applied. */
    std::vector<T> equilibrate() {
        std::vector<T> scale(n, num_traits<T>::from_int(1));
        for (std::size_t i = 0; i < n; ++i) {
            T m = num_traits<T>::from_int(0);
            for (std::size_t p = rowPtr[i]; p < rowPtr[i + 1]; ++p) {
                const T a = num_abs(T(val[p]));
                if (a > m) m = a;
            }
            if (m == num_traits<T>::from_int(0)) continue;
            scale[i] = m;
            if (m == num_traits<T>::from_int(1)) continue;
            for (std::size_t p = rowPtr[i]; p < rowPtr[i + 1]; ++p) val[p] /= m;
        }
        return scale;
    }

    /** Reordering so that entry (i,j) becomes (iperm[i], iperm[j]). */
    CsrMatrix permute_symmetric(const std::vector<std::size_t>& perm,
                                const std::vector<std::size_t>& iperm) const {
        CsrMatrix c;
        c.n = n;
        c.rowPtr.assign(n + 1, 0);
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t oi = perm[i];
            c.rowPtr[i + 1] = c.rowPtr[i] + (rowPtr[oi + 1] - rowPtr[oi]);
        }
        c.colIdx.resize(rowPtr[n]);
        c.val.resize(rowPtr[n]);
        std::vector<std::pair<std::size_t, T>> row;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t oi = perm[i];
            row.clear();
            for (std::size_t p = rowPtr[oi]; p < rowPtr[oi + 1]; ++p)
                row.push_back(std::make_pair(iperm[colIdx[p]], val[p]));
            std::sort(row.begin(), row.end(),
                      [](const std::pair<std::size_t, T>& a, const std::pair<std::size_t, T>& b) {
                          return a.first < b.first;
                      });
            std::size_t base = c.rowPtr[i];
            for (std::size_t k = 0; k < row.size(); ++k) {
                c.colIdx[base + k] = row[k].first;
                c.val[base + k] = row[k].second;
            }
        }
        c.build_diag();
        return c;
    }
};

/**
 * Reverse Cuthill-McKee ordering of the symmetrized pattern. Level structures
 * are grown from the lowest-degree unvisited node of each component, newly
 * discovered neighbours are appended in order of increasing degree, and the
 * result is reversed.
 */
template <class T>
std::vector<std::size_t> rcm_order(const CsrMatrix<T>& a) {
    const std::size_t n = a.n;
    std::vector<std::size_t> deg(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t p = a.rowPtr[i]; p < a.rowPtr[i + 1]; ++p) {
            const std::size_t j = a.colIdx[p];
            if (j == i) continue;
            ++deg[i];
            ++deg[j];
        }
    std::vector<std::size_t> adjPtr(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i) adjPtr[i + 1] = adjPtr[i] + deg[i];
    std::vector<std::size_t> adj(adjPtr[n]), fill(adjPtr.begin(), adjPtr.begin() + n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t p = a.rowPtr[i]; p < a.rowPtr[i + 1]; ++p) {
            const std::size_t j = a.colIdx[p];
            if (j == i) continue;
            adj[fill[i]++] = j;
            adj[fill[j]++] = i;
        }
    // Duplicate neighbours, from a pattern with both (i,j) and (j,i), only bias
    // the degree used for tie-breaking, which is harmless.

    std::vector<char> seen(n, 0);
    std::vector<std::size_t> result(n), queue(n);
    std::size_t count = 0;
    while (count < n) {
        std::size_t start = n;
        for (std::size_t i = 0; i < n; ++i)
            if (!seen[i] && (start == n || deg[i] < deg[start])) start = i;
        std::size_t head = count, tail = count;
        queue[tail++] = start;
        seen[start] = 1;
        while (head < tail) {
            const std::size_t v = queue[head++];
            result[count++] = v;
            const std::size_t from = tail;
            for (std::size_t p = adjPtr[v]; p < adjPtr[v + 1]; ++p) {
                const std::size_t u = adj[p];
                if (!seen[u]) {
                    seen[u] = 1;
                    queue[tail++] = u;
                }
            }
            std::stable_sort(queue.begin() + from, queue.begin() + tail,
                             [&deg](std::size_t x, std::size_t y) { return deg[x] < deg[y]; });
        }
    }
    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = result[n - 1 - i];
    return perm;
}

/** Threshold incomplete LU factorization, ILUT(p, tau) of Saad (1994). */
template <class T>
struct Ilut {
    bool ok = false;
    std::size_t n = 0;
    std::vector<std::size_t> lPtr, lCol, uPtr, uCol;
    std::vector<T> lVal, uVal, dInv;

    /** out = U^-1 L^-1 v, with L unit lower triangular. */
    void solve(const std::vector<T>& v, std::vector<T>& out) const {
        for (std::size_t i = 0; i < n; ++i) {
            T s = v[i];
            for (std::size_t p = lPtr[i]; p < lPtr[i + 1]; ++p) s -= lVal[p] * out[lCol[p]];
            out[i] = s;
        }
        for (std::size_t i = n; i-- > 0;) {
            T s = out[i];
            for (std::size_t p = uPtr[i]; p < uPtr[i + 1]; ++p) s -= uVal[p] * out[uCol[p]];
            out[i] = s * dInv[i];
        }
    }

    /** Keeps the entries of largest magnitude, by partial selection. */
    static std::size_t keep_largest(std::vector<std::size_t>& cols, const std::vector<T>& w,
                                    std::size_t len, std::size_t keep) {
        if (len <= keep) return len;
        for (std::size_t i = 0; i < keep; ++i) {
            std::size_t best = i;
            for (std::size_t j = i + 1; j < len; ++j)
                if (num_abs(T(w[cols[j]])) > num_abs(T(w[cols[best]]))) best = j;
            std::swap(cols[i], cols[best]);
        }
        return keep;
    }

    static Ilut factorize(const CsrMatrix<T>& a, double dropTol, double fillFactor) {
        Ilut f;
        const std::size_t n = a.n;
        f.n = n;
        const T zero = num_traits<T>::from_int(0);
        const std::size_t nnz = a.rowPtr[n];
        std::size_t lfil = static_cast<std::size_t>(
            std::ceil(fillFactor * static_cast<double>(nnz) / static_cast<double>(n == 0 ? 1 : n)));
        if (lfil < 1) lfil = 1;
        const T dropT = num_traits<T>::from_double(dropTol);

        f.lPtr.assign(n + 1, 0);
        f.uPtr.assign(n + 1, 0);
        f.dInv.assign(n, zero);

        std::vector<T> w(n, zero);
        std::vector<long> wPos(n, -1);
        std::vector<std::size_t> wIdx(n), rowsL(n), rowsU(n);
        std::size_t wCount = 0;
        std::priority_queue<std::size_t, std::vector<std::size_t>, std::greater<std::size_t>> pending;

        for (std::size_t i = 0; i < n; ++i) {
            T tnorm = zero;
            const std::size_t rowLen = a.rowPtr[i + 1] - a.rowPtr[i];
            for (std::size_t p = a.rowPtr[i]; p < a.rowPtr[i + 1]; ++p) {
                const std::size_t j = a.colIdx[p];
                w[j] = a.val[p];
                wPos[j] = static_cast<long>(wCount);
                wIdx[wCount++] = j;
                tnorm += num_abs(T(a.val[p]));
                if (j < i) pending.push(j);
            }
            if (rowLen == 0 || tnorm == zero) {
                // empty-row breakdown rationale: see _kb/03-api-layer.md (cpp port notes: mc)
                return f;  // ok stays false
            }
            const T tau = dropT * tnorm / num_traits<T>::from_int(static_cast<long>(rowLen));

            while (!pending.empty()) {
                const std::size_t k = pending.top();
                pending.pop();
                const T mult = w[k] * f.dInv[k];
                if (!(num_abs(T(mult)) > tau)) {
                    w[k] = zero;
                    continue;
                }
                w[k] = mult;
                for (std::size_t p = f.uPtr[k]; p < f.uPtr[k + 1]; ++p) {
                    const std::size_t j = f.uCol[p];
                    const T upd = mult * f.uVal[p];
                    if (wPos[j] >= 0) {
                        w[j] -= upd;
                    } else {
                        if (!(num_abs(T(upd)) > tau)) continue;
                        w[j] = -upd;
                        wPos[j] = static_cast<long>(wCount);
                        wIdx[wCount++] = j;
                        if (j < i) pending.push(j);
                    }
                }
            }

            std::size_t nl = 0, nu = 0;
            T diag = w[i];
            for (std::size_t t = 0; t < wCount; ++t) {
                const std::size_t j = wIdx[t];
                if (j == i) continue;
                if (!(num_abs(T(w[j])) > tau)) continue;
                if (j < i)
                    rowsL[nl++] = j;
                else
                    rowsU[nu++] = j;
            }
            nl = keep_largest(rowsL, w, nl, lfil);
            nu = keep_largest(rowsU, w, nu, lfil);
            std::sort(rowsL.begin(), rowsL.begin() + nl);
            std::sort(rowsU.begin(), rowsU.begin() + nu);

            for (std::size_t t = 0; t < nl; ++t) {
                f.lCol.push_back(rowsL[t]);
                f.lVal.push_back(w[rowsL[t]]);
            }
            for (std::size_t t = 0; t < nu; ++t) {
                f.uCol.push_back(rowsU[t]);
                f.uVal.push_back(w[rowsU[t]]);
            }
            f.lPtr[i + 1] = f.lCol.size();
            f.uPtr[i + 1] = f.uCol.size();

            // Saad's pivot remedy: see _kb/03-api-layer.md (cpp port notes: mc)
            if (!(num_abs(T(diag)) > tau) || !(diag == diag)) {
                const T substitute = tau > zero ? tau : num_traits<T>::from_double(1e-8);
                diag = diag < zero ? T(-substitute) : substitute;
            }
            f.dInv[i] = num_traits<T>::from_int(1) / diag;
            if (!num_isfinite(f.dInv[i])) return f;  // ok stays false

            for (std::size_t t = 0; t < wCount; ++t) {
                w[wIdx[t]] = zero;
                wPos[wIdx[t]] = -1;
            }
            wCount = 0;
            while (!pending.empty()) pending.pop();
        }
        f.ok = true;
        return f;
    }
};

/** ILUT if it factorized, Jacobi otherwise. */
template <class T>
struct Precond {
    Ilut<T> lu;
    std::vector<T> dinv;  ///< non-empty when the Jacobi fallback is in use

    static Precond of(const CsrMatrix<T>& a) {
        Precond m;
        m.lu = Ilut<T>::factorize(a, GMRES_ILUT_DROP_TOL, GMRES_ILUT_FILL_FACTOR);
        if (m.lu.ok) return m;
        const T zero = num_traits<T>::from_int(0);
        m.dinv.assign(a.n, num_traits<T>::from_int(1));
        for (std::size_t i = 0; i < a.n; ++i) {
            const T d = a.diagPtr[i] < 0 ? zero : a.val[static_cast<std::size_t>(a.diagPtr[i])];
            m.dinv[i] = (d == zero) ? num_traits<T>::from_int(1) : T(num_traits<T>::from_int(1) / d);
        }
        return m;
    }

    void apply(const std::vector<T>& v, std::vector<T>& out) const {
        if (!dinv.empty()) {
            for (std::size_t i = 0; i < dinv.size(); ++i) out[i] = dinv[i] * v[i];
            return;
        }
        lu.solve(v, out);
    }
};

/** Equilibrated, reordered and preconditioned form of a coefficient matrix. */
template <class T>
struct GmresPrepared {
    std::size_t n = 0;
    CsrMatrix<T> csr;
    std::vector<std::size_t> perm;
    std::vector<T> rowScale;
    Precond<T> M;

    explicit GmresPrepared(const Matrix<T>& A) {
        n = A.rows();
        if (A.cols() != n) throw InputError("ctmc_gmres: matrix is not square");
        CsrMatrix<T> c = CsrMatrix<T>::of(A);
        rowScale = c.equilibrate();
        perm = rcm_order(c);
        std::vector<std::size_t> iperm(n);
        for (std::size_t i = 0; i < n; ++i) iperm[perm[i]] = i;
        csr = c.permute_symmetric(perm, iperm);
        M = Precond<T>::of(csr);
    }
};

template <class T>
T vec_norm2(const std::vector<T>& v) {
    T s = num_traits<T>::from_int(0);
    for (const T& x : v) s += x * x;
    using std::sqrt;
    return sqrt(s);
}

template <class T>
T vec_dot(const std::vector<T>& a, const std::vector<T>& b) {
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

/** Right-preconditioned restarted GMRES on an already prepared system. */
template <class T>
GmresResult<T> gmres_solve(const GmresPrepared<T>& prep, const std::vector<T>& rhsIn,
                           const std::vector<T>& x0In, double tol, long restart, long maxit) {
    const std::size_t n = prep.n;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (tol <= 0.0) tol = GMRES_DEFAULT_TOL;
    if (restart <= 0) restart = std::min(static_cast<long>(n), GMRES_DEFAULT_RESTART);
    restart = std::min(restart, static_cast<long>(n));
    if (maxit <= 0)
        maxit = static_cast<long>(std::ceil(static_cast<double>(n) / static_cast<double>(restart)));
    maxit = std::max(1L, std::min(maxit, static_cast<long>(n)));
    const std::size_t m = static_cast<std::size_t>(restart);
    const T tolT = num_traits<T>::from_double(tol);

    // Row scaling, then RCM permutation, on both the right-hand side and the
    // initial guess.
    std::vector<T> rhs(n), x(n);
    for (std::size_t i = 0; i < n; ++i) rhs[i] = rhsIn[prep.perm[i]] / prep.rowScale[prep.perm[i]];
    for (std::size_t i = 0; i < n; ++i) x[i] = x0In[prep.perm[i]];

    T bnorm = vec_norm2(rhs);
    if (bnorm == zero) bnorm = one;

    std::vector<T> r(n), w(n), z(n), corr(n);
    prep.csr.mult(x, r);
    for (std::size_t i = 0; i < n; ++i) r[i] = rhs[i] - r[i];
    T beta = vec_norm2(r);

    GmresResult<T> out;
    out.flag = 1;
    out.iter = 0;
    out.relres = beta / bnorm;
    if (out.relres <= tolT) {
        out.flag = 0;
        out.x.assign(n, zero);
        for (std::size_t i = 0; i < n; ++i) out.x[prep.perm[i]] = x[i];
        return out;
    }

    std::vector<std::vector<T>> V(m + 1, std::vector<T>(n, zero));
    std::vector<std::vector<T>> H(m + 1, std::vector<T>(m, zero));
    std::vector<T> cs(m, zero), sn(m, zero), g(m + 1, zero), y(m, zero);
    const T breakTol = num_traits<T>::from_double(GMRES_BREAKDOWN_TOL);
    const T initres = out.relres;
    bool havePrev = false;
    T prevrelres = zero;

    for (long cycle = 0; cycle < maxit; ++cycle) {
        beta = vec_norm2(r);
        if (beta == zero) {
            out.flag = 0;
            out.relres = zero;
            break;
        }
        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i] / beta;
        std::fill(g.begin(), g.end(), zero);
        g[0] = beta;

        std::size_t k = 0;
        for (std::size_t j = 0; j < m; ++j) {
            // right-preconditioning rationale: see _kb/03-api-layer.md (cpp port notes: mc)
            prep.M.apply(V[j], z);
            prep.csr.mult(z, w);
            ++out.iter;

            const T wnorm0 = vec_norm2(w);
            // Gram-Schmidt reorthogonalization rationale: see _kb/03-api-layer.md (cpp port notes: mc)
            for (int pass = 0; pass < 2; ++pass)
                for (std::size_t i = 0; i <= j; ++i) {
                    const T hij = vec_dot(V[i], w);
                    H[i][j] += hij;
                    for (std::size_t q = 0; q < n; ++q) w[q] -= hij * V[i][q];
                }
            const T hnext = vec_norm2(w);
            H[j + 1][j] = hnext;

            k = j + 1;
            const bool breakdown = !(hnext > breakTol * wnorm0);
            if (!breakdown)
                for (std::size_t q = 0; q < n; ++q) V[j + 1][q] = w[q] / hnext;

            // Accumulated Givens rotations on the new Hessenberg column, then a
            // fresh rotation annihilating its subdiagonal entry.
            for (std::size_t i = 0; i < j; ++i) {
                const T t1 = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                H[i][j] = t1;
            }
            using std::sqrt;
            const T denom = sqrt(T(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]));
            if (denom == zero) {
                cs[j] = one;
                sn[j] = zero;
            } else {
                cs[j] = H[j][j] / denom;
                sn[j] = H[j + 1][j] / denom;
            }
            H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
            H[j + 1][j] = zero;
            g[j + 1] = -sn[j] * g[j];
            g[j] = cs[j] * g[j];

            out.relres = num_abs(T(g[j + 1])) / bnorm;
            if (out.relres <= tolT || breakdown) break;
        }

        // Least-squares solution on the rotated Hessenberg system, mapped back
        // through the preconditioner.
        for (std::size_t i = k; i-- > 0;) {
            T s = g[i];
            for (std::size_t q = i + 1; q < k; ++q) s -= H[i][q] * y[q];
            y[i] = (H[i][i] == zero) ? zero : T(s / H[i][i]);
        }
        std::fill(corr.begin(), corr.end(), zero);
        for (std::size_t i = 0; i < k; ++i)
            for (std::size_t q = 0; q < n; ++q) corr[q] += y[i] * V[i][q];
        prep.M.apply(corr, z);
        for (std::size_t q = 0; q < n; ++q) x[q] += z[q];

        prep.csr.mult(x, r);
        for (std::size_t q = 0; q < n; ++q) r[q] = rhs[q] - r[q];
        out.relres = vec_norm2(r) / bnorm;

        for (std::size_t i = 0; i <= m; ++i) std::fill(H[i].begin(), H[i].end(), zero);

        if (out.relres <= tolT) {
            out.flag = 0;
            break;
        }
        // divergence-detection rationale: see _kb/03-api-layer.md (cpp port notes: mc)
        if (!(out.relres < num_traits<T>::from_int(100) * initres)) {
            out.flag = 3;
            break;
        }
        // Stagnation: a cycle that fails to reduce the residual will not do so
        // on the next one either.
        if (havePrev && out.relres >= prevrelres * (one - num_traits<T>::from_double(1e-12))) {
            out.flag = 3;
            break;
        }
        prevrelres = out.relres;
        havePrev = true;
    }

    out.x.assign(n, zero);
    for (std::size_t i = 0; i < n; ++i) out.x[prep.perm[i]] = x[i];
    for (std::size_t i = 0; i < n; ++i)
        if (!num_isfinite(out.x[i])) {
            out.flag = 3;
            out.relres = num_traits<T>::from_double(1e300);
            return out;
        }
    if (out.relres <= tolT) out.flag = 0;
    return out;
}

}  // namespace detail

/**
 * @param A       coefficient matrix, already assembled
 * @param b       right-hand side
 * @param tol     relative residual tolerance (default 1e-12)
 * @param restart restart length; <= 0 selects min(n, 50)
 * @param maxit   outer cycles; <= 0 selects ceil(n / restart)
 * @param x0      initial guess; empty selects the uniform vector ones(n)/n
 */
template <class T>
GmresResult<T> ctmc_gmres(const Matrix<T>& A, const std::vector<T>& b, double tol = 1e-12,
                          long restart = 0, long maxit = 0, const std::vector<T>& x0 = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_gmres requires transcendental arithmetic: the Arnoldi step normalizes by a "
                  "Euclidean norm and the iteration stops on a residual tolerance, so there is no "
                  "exact result to converge to; use ctmc_solve for an exact solve");
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("ctmc_gmres: matrix is not square");
    if (b.size() != n) throw InputError("ctmc_gmres: right-hand side has the wrong length");
    if (!x0.empty() && x0.size() != n) throw InputError("ctmc_gmres: initial guess has the wrong length");
    std::vector<T> guess = x0;
    if (guess.empty())
        guess.assign(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(n)));
    const detail::GmresPrepared<T> prep(A);
    return detail::gmres_solve(prep, b, guess, tol, restart, maxit);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_GMRES_H
