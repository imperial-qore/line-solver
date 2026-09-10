/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_SIMPLEX_H
#define LINE_UTIL_SIMPLEX_H

/**
 * Templated primal simplex with Bland's rule.
 *
 * Solves
 *      max (or min)  c'x
 *      subject to    A x <= b,  Aeq x = beq,  l <= x <= u
 * where l and u are per-variable and either bound may be absent.
 *
 * Why this exists: the mapqn quadratic-reduction bounds are linear programs
 * whose optimum IS the bound being reported. MATLAB reaches it with
 * interior-point linprog and lands a few digits short (see the MATLAB accuracy
 * note in matlab/lib/qrf/mapqn_bnd_qr_ld.m); the JAR reaches it with Apache
 * Commons SimplexSolver in double precision. Instantiated at line::Rational
 * this solver is EXACT: rational data implies a rational optimum, every pivot
 * is a field operation on rationals, and the reported bound is the vertex
 * value with no rounding anywhere. That is the point of the port, so the
 * pivoting rule must not depend on a tolerance.
 *
 * Bland's rule is what makes that possible. It selects the lowest-index
 * column with a strictly favourable reduced cost and, among the rows attaining
 * the minimum ratio, the one whose basic variable has the lowest index. That
 * pair of rules is enough to prove finite termination with no anti-cycling
 * perturbation and no tolerance, which is why the solver is NOT gated on
 * num_traits<T>::has_transcendental: it uses only +, -, *, / and comparison,
 * all of which Rational supports exactly. At inexact T a small tolerance is
 * used for the sign tests, purely so that round-off does not report a
 * favourable reduced cost that is really zero; at exact T the tolerance is
 * exactly zero and no such fudge exists.
 *
 * Variable bounds are handled by the solver itself, not by the caller:
 *   - l_j == u_j            the variable is substituted out (fixed), which is
 *                           what makes the mapqn assembly tractable, since its
 *                           ZERO1/2/3 families fix the majority of variables
 *                           at zero;
 *   - l_j finite            x_j = l_j + y_j with y_j >= 0, and a finite u_j
 *                           becomes one extra row y_j <= u_j - l_j;
 *   - l_j absent, u_j given x_j = u_j - y_j with y_j >= 0;
 *   - both absent           x_j = y_j^+ - y_j^- with both parts >= 0.
 * Callers therefore never need to add explicit 0 <= x <= 1 rows the way the
 * JAR must for Apache SimplexSolver (see the mapqn note in _kb/03-api-layer.md);
 * they call set_bounds and the rows appear internally only where a finite
 * upper bound actually needs one.
 *
 * Assembly is sparse. LpModel accumulates a row in a dense scratch vector with
 * a touched-index list and emits only the nonzeros, so building the mapqn LPs
 * costs O(nnz) memory rather than the O(rows*cols) a dense builder would need.
 * Repeated add() on the same column accumulates, matching the
 * `row(idx) = row(idx) + v` idiom the MATLAB reference uses.
 *
 * The tableau itself is dense: it holds B^{-1}[A I] explicitly. This is a
 * deliberate scope choice -- the port targets small and medium instances where
 * exactness is the objective, not the large blocking instances that need a
 * sparse revised simplex with LU updates.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace lp {

/** Outcome of a solve. */
enum class LpStatus {
    Optimal,        ///< an optimal vertex was reached
    Infeasible,     ///< phase 1 ended with residual artificial mass
    Unbounded,      ///< an improving column has no blocking row
    IterationLimit  ///< the iteration cap was hit (cannot happen under Bland's rule
                    ///< with exact arithmetic; a guard for inexact T)
};

inline const char* lp_status_name(LpStatus s) {
    switch (s) {
        case LpStatus::Optimal: return "Optimal";
        case LpStatus::Infeasible: return "Infeasible";
        case LpStatus::Unbounded: return "Unbounded";
        default: return "IterationLimit";
    }
}

template <class T>
struct LpSolution {
    LpStatus status = LpStatus::Infeasible;
    std::vector<T> x;   ///< primal solution in the ORIGINAL variable space
    T objective = T();  ///< c'x, in the sense requested (max or min)
    std::size_t iterations = 0;
    bool ok() const { return status == LpStatus::Optimal; }
};

/** Row relation. */
enum class LpSense { LE, EQ, GE };

/**
 * Sparse LP in the natural form, with per-variable bounds.
 *
 * Default bounds are x_j >= 0 with no upper bound, matching linprog's
 * convention when lb is given as zeros and ub is omitted.
 */
template <class T>
class LpModel {
public:
    explicit LpModel(std::size_t nvars)
        : n_(nvars),
          lb_(nvars, T()),
          ub_(nvars, T()),
          lb_inf_(nvars, 0),
          ub_inf_(nvars, 1),
          c_(nvars, T()),
          scratch_(nvars, T()),
          touched_flag_(nvars, 0) {}

    std::size_t num_vars() const { return n_; }
    std::size_t num_rows() const { return rhs_.size(); }
    std::size_t num_nonzeros() const { return cols_.size(); }

    // ---------------------------------------------------------------- bounds
    void set_lower(std::size_t j, const T& v) {
        check(j);
        lb_[j] = v;
        lb_inf_[j] = 0;
    }
    void set_upper(std::size_t j, const T& v) {
        check(j);
        ub_[j] = v;
        ub_inf_[j] = 0;
    }
    void set_bounds(std::size_t j, const T& lo, const T& hi) {
        set_lower(j, lo);
        set_upper(j, hi);
    }
    void set_free_lower(std::size_t j) {
        check(j);
        lb_inf_[j] = 1;
    }
    void set_free_upper(std::size_t j) {
        check(j);
        ub_inf_[j] = 1;
    }
    void set_free(std::size_t j) {
        set_free_lower(j);
        set_free_upper(j);
    }
    /** Pin a variable to a value; it is substituted out of the tableau. */
    void fix(std::size_t j, const T& v) { set_bounds(j, v, v); }

    const T& lower(std::size_t j) const { return lb_[j]; }
    const T& upper(std::size_t j) const { return ub_[j]; }
    bool lower_is_free(std::size_t j) const { return lb_inf_[j] != 0; }
    bool upper_is_free(std::size_t j) const { return ub_inf_[j] != 0; }

    // ------------------------------------------------------------- objective
    void set_cost(std::size_t j, const T& v) {
        check(j);
        c_[j] = v;
    }
    void add_cost(std::size_t j, const T& v) {
        check(j);
        c_[j] += v;
    }
    const std::vector<T>& costs() const { return c_; }
    /** true to maximize c'x (the default), false to minimize. */
    void set_maximize(bool m) { maximize_ = m; }
    bool maximize() const { return maximize_; }

    // ---------------------------------------------------------- row assembly
    /** Discard whatever the row accumulator holds. */
    void row_clear() {
        for (std::size_t k = 0; k < touched_.size(); ++k) {
            scratch_[touched_[k]] = T();
            touched_flag_[touched_[k]] = 0;
        }
        touched_.clear();
    }

    /** row(j) += v, the accumulation the MATLAB reference performs. */
    void row_add(std::size_t j, const T& v) {
        check(j);
        if (!touched_flag_[j]) {
            touched_flag_[j] = 1;
            touched_.push_back(j);
        }
        scratch_[j] += v;
    }

    void row_add_int(std::size_t j, long v) { row_add(j, num_traits<T>::from_int(v)); }

    /** Emit the accumulated row with the given relation and right-hand side. */
    void emit(LpSense sense, const T& rhs) {
        const T zero = T();
        std::size_t nnz = 0;
        for (std::size_t k = 0; k < touched_.size(); ++k) {
            const std::size_t j = touched_[k];
            if (scratch_[j] != zero) {
                cols_.push_back(j);
                vals_.push_back(scratch_[j]);
                ++nnz;
            }
        }
        row_start_.push_back(row_start_.back() + nnz);
        rhs_.push_back(rhs);
        sense_.push_back(sense);
        row_clear();
    }
    void emit_le(const T& rhs) { emit(LpSense::LE, rhs); }
    void emit_eq(const T& rhs) { emit(LpSense::EQ, rhs); }
    void emit_ge(const T& rhs) { emit(LpSense::GE, rhs); }
    void emit_le_int(long rhs) { emit(LpSense::LE, num_traits<T>::from_int(rhs)); }
    void emit_eq_int(long rhs) { emit(LpSense::EQ, num_traits<T>::from_int(rhs)); }
    void emit_ge_int(long rhs) { emit(LpSense::GE, num_traits<T>::from_int(rhs)); }

    // --------------------------------------------------------- row inspection
    std::size_t row_begin(std::size_t i) const { return row_start_[i]; }
    std::size_t row_end(std::size_t i) const { return row_start_[i + 1]; }
    std::size_t col_at(std::size_t k) const { return cols_[k]; }
    const T& val_at(std::size_t k) const { return vals_[k]; }
    const T& rhs(std::size_t i) const { return rhs_[i]; }
    LpSense sense(std::size_t i) const { return sense_[i]; }

private:
    void check(std::size_t j) const {
        if (j >= n_) throw InputError("LpModel: variable index out of range");
    }

    std::size_t n_;
    std::vector<T> lb_, ub_;
    std::vector<char> lb_inf_, ub_inf_;
    std::vector<T> c_;
    bool maximize_ = true;

    std::vector<std::size_t> row_start_ = std::vector<std::size_t>(1, 0);
    std::vector<std::size_t> cols_;
    std::vector<T> vals_;
    std::vector<T> rhs_;
    std::vector<LpSense> sense_;

    std::vector<T> scratch_;
    std::vector<std::size_t> touched_;
    std::vector<char> touched_flag_;
};

// ---------------------------------------------------------------------------
// Sign tests. Exactly zero tolerance when T is exact, which is the whole point.
// ---------------------------------------------------------------------------

template <class T>
inline T simplex_tolerance() {
    return num_traits<T>::is_exact ? T() : num_traits<T>::from_double(1e-10);
}

namespace detail {

/** Internal standard-form column: how an original variable maps into y >= 0. */
enum class VarKind {
    Fixed,      ///< l == u, substituted out
    ShiftUp,    ///< x = l + y
    ShiftDown,  ///< x = u - y   (no lower bound, finite upper bound)
    Free        ///< x = y+ - y-
};

struct VarMap {
    VarKind kind = VarKind::ShiftUp;
    std::size_t pos = 0;  ///< column of y (or y+)
    std::size_t neg = 0;  ///< column of y- when Free
};

}  // namespace detail

/**
 * Solve the model. See the header comment for the bound handling; the returned
 * x is in the caller's variable space and satisfies the bounds exactly when T
 * is exact.
 */
template <class T>
LpSolution<T> simplex_solve(const LpModel<T>& model, std::size_t max_iterations = 0) {
    using detail::VarKind;
    using detail::VarMap;

    const T zero = T();
    const T one = num_traits<T>::from_int(1);
    const T tol = simplex_tolerance<T>();
    const std::size_t n = model.num_vars();

    // ---- variable mapping -------------------------------------------------
    std::vector<VarMap> vmap(n);
    std::vector<T> fixed_value(n, zero);
    std::vector<T> offset(n, zero);  // the constant part of x_j
    std::size_t ny = 0;
    for (std::size_t j = 0; j < n; ++j) {
        const bool lf = model.lower_is_free(j), uf = model.upper_is_free(j);
        if (!lf && !uf && model.lower(j) == model.upper(j)) {
            vmap[j].kind = VarKind::Fixed;
            fixed_value[j] = model.lower(j);
            offset[j] = model.lower(j);
        } else if (!lf && !uf && model.upper(j) < model.lower(j)) {
            LpSolution<T> s;
            s.status = LpStatus::Infeasible;
            return s;
        } else if (!lf) {
            vmap[j].kind = VarKind::ShiftUp;
            vmap[j].pos = ny++;
            offset[j] = model.lower(j);
        } else if (!uf) {
            vmap[j].kind = VarKind::ShiftDown;
            vmap[j].pos = ny++;
            offset[j] = model.upper(j);
        } else {
            vmap[j].kind = VarKind::Free;
            vmap[j].pos = ny++;
            vmap[j].neg = ny++;
            offset[j] = zero;
        }
    }

    // ---- standard-form rows: G y <= g and H y = h --------------------------
    // Each stored as (cols, vals, rhs). Inequalities first, then equalities.
    std::vector<std::vector<std::size_t>> le_cols, eq_cols;
    std::vector<std::vector<T>> le_vals, eq_vals;
    std::vector<T> le_rhs, eq_rhs;

    for (std::size_t i = 0; i < model.num_rows(); ++i) {
        const LpSense sn = model.sense(i);
        const bool negate = (sn == LpSense::GE);
        std::vector<std::size_t> cs;
        std::vector<T> vs;
        T r = model.rhs(i);
        for (std::size_t k = model.row_begin(i); k < model.row_end(i); ++k) {
            const std::size_t j = model.col_at(k);
            const T a = model.val_at(k);
            if (vmap[j].kind == VarKind::Fixed) {
                const T contrib = a * fixed_value[j];
                r -= contrib;
                continue;
            }
            if (offset[j] != zero) {
                const T contrib = a * offset[j];
                r -= contrib;
            }
            if (vmap[j].kind == VarKind::ShiftUp) {
                cs.push_back(vmap[j].pos);
                vs.push_back(a);
            } else if (vmap[j].kind == VarKind::ShiftDown) {
                cs.push_back(vmap[j].pos);
                const T na = -a;
                vs.push_back(na);
            } else {
                cs.push_back(vmap[j].pos);
                vs.push_back(a);
                cs.push_back(vmap[j].neg);
                const T na = -a;
                vs.push_back(na);
            }
        }
        if (negate) {
            for (std::size_t k = 0; k < vs.size(); ++k) {
                const T nv = -vs[k];
                vs[k] = nv;
            }
            const T nr = -r;
            r = nr;
        }
        if (sn == LpSense::EQ) {
            eq_cols.push_back(cs);
            eq_vals.push_back(vs);
            eq_rhs.push_back(r);
        } else {
            le_cols.push_back(cs);
            le_vals.push_back(vs);
            le_rhs.push_back(r);
        }
    }

    // Finite upper bounds on shifted variables become one row each.
    for (std::size_t j = 0; j < n; ++j) {
        if (vmap[j].kind == VarKind::ShiftUp && !model.upper_is_free(j)) {
            const T span = model.upper(j) - model.lower(j);
            le_cols.push_back(std::vector<std::size_t>(1, vmap[j].pos));
            le_vals.push_back(std::vector<T>(1, one));
            le_rhs.push_back(span);
        } else if (vmap[j].kind == VarKind::ShiftDown && !model.lower_is_free(j)) {
            const T span = model.upper(j) - model.lower(j);
            le_cols.push_back(std::vector<std::size_t>(1, vmap[j].pos));
            le_vals.push_back(std::vector<T>(1, one));
            le_rhs.push_back(span);
        }
    }

    const std::size_t n_le = le_rhs.size();
    const std::size_t n_eq = eq_rhs.size();
    std::size_t m = n_le + n_eq;
    const std::size_t n_struct = ny + n_le;  // structural + slack columns

    // ---- objective in y space ---------------------------------------------
    std::vector<T> cy(n_struct, zero);
    T const_obj = zero;
    for (std::size_t j = 0; j < n; ++j) {
        const T cj = model.costs()[j];
        if (cj == zero) continue;
        if (vmap[j].kind == VarKind::Fixed) {
            const T contrib = cj * fixed_value[j];
            const_obj += contrib;
            continue;
        }
        if (offset[j] != zero) {
            const T contrib = cj * offset[j];
            const_obj += contrib;
        }
        if (vmap[j].kind == VarKind::ShiftUp) {
            cy[vmap[j].pos] += cj;
        } else if (vmap[j].kind == VarKind::ShiftDown) {
            cy[vmap[j].pos] -= cj;
        } else {
            cy[vmap[j].pos] += cj;
            cy[vmap[j].neg] -= cj;
        }
    }
    if (!model.maximize()) {
        for (std::size_t k = 0; k < cy.size(); ++k) {
            const T nv = -cy[k];
            cy[k] = nv;
        }
    }

    // tableau layout ([y | slacks | artificials]): see _kb/14-cpp-multiprecision.md
    const std::size_t ncol1 = n_struct + m;
    std::vector<T> tab(m * (ncol1 + 1), zero);
    const std::size_t stride1 = ncol1 + 1;
    for (std::size_t i = 0; i < n_le; ++i) {
        for (std::size_t k = 0; k < le_cols[i].size(); ++k)
            tab[i * stride1 + le_cols[i][k]] += le_vals[i][k];
        tab[i * stride1 + ny + i] = one;  // slack
        tab[i * stride1 + ncol1] = le_rhs[i];
    }
    for (std::size_t i = 0; i < n_eq; ++i) {
        const std::size_t r = n_le + i;
        for (std::size_t k = 0; k < eq_cols[i].size(); ++k)
            tab[r * stride1 + eq_cols[i][k]] += eq_vals[i][k];
        tab[r * stride1 + ncol1] = eq_rhs[i];
    }
    // Nonnegative right-hand sides, then an artificial basis.
    for (std::size_t i = 0; i < m; ++i) {
        if (tab[i * stride1 + ncol1] < zero) {
            for (std::size_t j = 0; j <= ncol1; ++j) {
                const T nv = -tab[i * stride1 + j];
                tab[i * stride1 + j] = nv;
            }
        }
        tab[i * stride1 + n_struct + i] = one;
    }

    std::vector<std::size_t> basis(m);
    for (std::size_t i = 0; i < m; ++i) basis[i] = n_struct + i;

    std::size_t iter_cap = max_iterations;
    if (iter_cap == 0) {
        const std::size_t base = (m + 1) * (ncol1 + 1);
        iter_cap = base < 100000 ? 100000 : base * 20;
    }
    std::size_t iters = 0;

    // ---- shared pivot loop (Bland's rule) ---------------------------------
    // cost points at a vector of length ncols; entering columns are restricted
    // to [0, ncols_allowed).
    struct Pivot {
        static void apply(std::vector<T>& tb, std::size_t stride, std::size_t rows, std::size_t row,
                          std::size_t col) {
            const T piv = tb[row * stride + col];
            const T inv = num_traits<T>::from_int(1) / piv;
            for (std::size_t j = 0; j < stride; ++j) {
                const T nv = tb[row * stride + j] * inv;
                tb[row * stride + j] = nv;
            }
            tb[row * stride + col] = num_traits<T>::from_int(1);
            for (std::size_t i = 0; i < rows; ++i) {
                if (i == row) continue;
                const T f = tb[i * stride + col];
                if (f == T()) continue;
                for (std::size_t j = 0; j < stride; ++j) {
                    const T nv = tb[i * stride + j] - f * tb[row * stride + j];
                    tb[i * stride + j] = nv;
                }
                tb[i * stride + col] = T();
            }
        }
    };

    // Phase 1: maximize -(sum of artificials).
    {
        std::vector<T> c1(ncol1, zero);
        for (std::size_t i = 0; i < m; ++i) c1[n_struct + i] = -one;
        bool unbounded = false;
        while (true) {
            if (++iters > iter_cap) {
                LpSolution<T> s;
                s.status = LpStatus::IterationLimit;
                return s;
            }
            // reduced costs d_j = c_j - cB' * col_j
            std::size_t enter = ncol1;
            for (std::size_t j = 0; j < ncol1; ++j) {
                T d = c1[j];
                for (std::size_t i = 0; i < m; ++i) {
                    const T cb = c1[basis[i]];
                    if (cb == zero) continue;
                    const T contrib = cb * tab[i * stride1 + j];
                    d -= contrib;
                }
                if (d > tol) {
                    enter = j;
                    break;  // Bland: lowest index
                }
            }
            if (enter == ncol1) break;  // phase-1 optimum
            std::size_t leave = m;
            T best_num = zero, best_den = one;
            for (std::size_t i = 0; i < m; ++i) {
                const T a = tab[i * stride1 + enter];
                if (!(a > tol)) continue;
                const T rr = tab[i * stride1 + ncol1];
                if (leave == m) {
                    leave = i;
                    best_num = rr;
                    best_den = a;
                } else {
                    const T lhs = rr * best_den;
                    const T rhs2 = best_num * a;
                    if (lhs < rhs2 || (lhs == rhs2 && basis[i] < basis[leave])) {
                        leave = i;
                        best_num = rr;
                        best_den = a;
                    }
                }
            }
            if (leave == m) {
                unbounded = true;  // cannot happen: phase-1 objective is bounded
                break;
            }
            Pivot::apply(tab, stride1, m, leave, enter);
            basis[leave] = enter;
        }
        if (unbounded) {
            LpSolution<T> s;
            s.status = LpStatus::Infeasible;
            return s;
        }
        // residual artificial mass -> infeasible
        T infeas = zero;
        for (std::size_t i = 0; i < m; ++i)
            if (basis[i] >= n_struct) infeas += tab[i * stride1 + ncol1];
        if (infeas > tol) {
            LpSolution<T> s;
            s.status = LpStatus::Infeasible;
            return s;
        }
        // Drive artificials out of the basis; rows with no pivot are redundant.
        std::vector<char> drop(m, 0);
        for (std::size_t i = 0; i < m; ++i) {
            if (basis[i] < n_struct) continue;
            std::size_t piv = n_struct;
            for (std::size_t j = 0; j < n_struct; ++j) {
                const T a = tab[i * stride1 + j];
                if (a > tol || a < -tol) {
                    piv = j;
                    break;
                }
            }
            if (piv == n_struct) {
                drop[i] = 1;
            } else {
                Pivot::apply(tab, stride1, m, i, piv);
                basis[i] = piv;
            }
        }
        // Rebuild the tableau without artificial columns and dropped rows.
        const std::size_t stride2 = n_struct + 1;
        std::vector<T> tab2;
        std::vector<std::size_t> basis2;
        tab2.reserve(m * stride2);
        for (std::size_t i = 0; i < m; ++i) {
            if (drop[i]) continue;
            for (std::size_t j = 0; j < n_struct; ++j) tab2.push_back(tab[i * stride1 + j]);
            tab2.push_back(tab[i * stride1 + ncol1]);
            basis2.push_back(basis[i]);
        }
        tab.swap(tab2);
        basis.swap(basis2);
        m = basis.size();
    }

    // Phase 2.
    const std::size_t stride = n_struct + 1;
    bool unbounded = false;
    while (true) {
        if (++iters > iter_cap) {
            LpSolution<T> s;
            s.status = LpStatus::IterationLimit;
            return s;
        }
        std::size_t enter = n_struct;
        for (std::size_t j = 0; j < n_struct; ++j) {
            T d = cy[j];
            for (std::size_t i = 0; i < m; ++i) {
                const T cb = cy[basis[i]];
                if (cb == zero) continue;
                const T contrib = cb * tab[i * stride + j];
                d -= contrib;
            }
            if (d > tol) {
                enter = j;
                break;
            }
        }
        if (enter == n_struct) break;
        std::size_t leave = m;
        T best_num = zero, best_den = one;
        for (std::size_t i = 0; i < m; ++i) {
            const T a = tab[i * stride + enter];
            if (!(a > tol)) continue;
            const T rr = tab[i * stride + n_struct];
            if (leave == m) {
                leave = i;
                best_num = rr;
                best_den = a;
            } else {
                const T lhs = rr * best_den;
                const T rhs2 = best_num * a;
                if (lhs < rhs2 || (lhs == rhs2 && basis[i] < basis[leave])) {
                    leave = i;
                    best_num = rr;
                    best_den = a;
                }
            }
        }
        if (leave == m) {
            unbounded = true;
            break;
        }
        Pivot::apply(tab, stride, m, leave, enter);
        basis[leave] = enter;
    }

    LpSolution<T> sol;
    sol.iterations = iters;
    if (unbounded) {
        sol.status = LpStatus::Unbounded;
        return sol;
    }
    sol.status = LpStatus::Optimal;

    std::vector<T> y(n_struct, zero);
    for (std::size_t i = 0; i < m; ++i)
        if (basis[i] < n_struct) y[basis[i]] = tab[i * stride + n_struct];

    sol.x.assign(n, zero);
    for (std::size_t j = 0; j < n; ++j) {
        switch (vmap[j].kind) {
            case VarKind::Fixed: sol.x[j] = fixed_value[j]; break;
            case VarKind::ShiftUp: sol.x[j] = offset[j] + y[vmap[j].pos]; break;
            case VarKind::ShiftDown: sol.x[j] = offset[j] - y[vmap[j].pos]; break;
            default: sol.x[j] = y[vmap[j].pos] - y[vmap[j].neg]; break;
        }
    }
    T obj = zero;
    for (std::size_t j = 0; j < n; ++j) {
        const T cj = model.costs()[j];
        if (cj == zero) continue;
        const T contrib = cj * sol.x[j];
        obj += contrib;
    }
    sol.objective = obj;
    (void)const_obj;
    return sol;
}

}  // namespace lp
}  // namespace line

#endif  // LINE_UTIL_SIMPLEX_H
