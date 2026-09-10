/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_LP_HIGHS_H
#define LINE_UTIL_LP_HIGHS_H

/**
 * A sparse LP backend for line::lp::LpModel, on HiGHS (MIT).
 *
 * WHY IT EXISTS. util/simplex.h carries a DENSE tableau. That is the right
 * choice for what it was written for -- it is exact under line::Rational, uses
 * Bland's rule with no tolerance, and needs no dependency -- but it costs
 * O(rows * cols) per pivot, so it clears a few hundred columns and stalls on a
 * few thousand. The QRF blocking bounds are MR * B^2 columns
 * (mapqn_qr_bounds_bas.h) and B^2 (mapqn_qr_bounds_rsrd.h), which is 3894 for
 * even example_bas_small.m and tens of thousands for the paper instances. This
 * backend is what makes those reachable.
 *
 * IT IS NOT A REPLACEMENT, AND MUST NOT BECOME ONE. HiGHS is double precision.
 * Four headers in this tree promise that at T = line::Rational the returned
 * value is the EXACT optimum of the exact polytope, and the mapqn tests assert
 * that as equalities on fractions (3/4, 2/3, 3/7, 6/7), not as tolerances.
 * Routing Rational here would silently turn those equalities into rounding.
 * So the dispatcher below refuses any T other than double, at compile time.
 *
 * PRESOLVE IS LEFT ON but the caller should know it exists: the QRF equality
 * blocks are heavily redundant (matlab/lib/qrf/qrf_independent_rows.m selects a
 * maximal independent subset by pivoted QR precisely because linprog struggles
 * otherwise), and presolve is what absorbs that redundancy here. If a model
 * ever comes back Infeasible where the dense solver says Optimal, re-run it
 * with presolve off before believing the answer -- that is the analogue of the
 * `adaptive_rho` trap recorded for the OSQP-backed bounds in _kb.
 */

#include <cstddef>
#include <vector>

#include "line/util/simplex.h"

#ifdef LINE_MP_HAVE_HIGHS
#include "Highs.h"
#endif

namespace line {
namespace lp {

/** True when a sparse backend is compiled in. */
inline bool highs_available() {
#ifdef LINE_MP_HAVE_HIGHS
    return true;
#else
    return false;
#endif
}

#ifdef LINE_MP_HAVE_HIGHS

/**
 * Solve a double-precision LpModel with HiGHS.
 *
 * The translation is direct: LpModel already stores rows in compressed sparse
 * row form, which is one of the two layouts HiGHS accepts, so the matrix is
 * handed over without a transpose. Free bounds become +-kHighsInf.
 */
inline LpSolution<double> highs_solve(const LpModel<double>& model) {
    LpSolution<double> out;
    const std::size_t n = model.num_vars(), m = model.num_rows();

    HighsModel hm;
    hm.lp_.num_col_ = static_cast<HighsInt>(n);
    hm.lp_.num_row_ = static_cast<HighsInt>(m);
    hm.lp_.sense_ = model.maximize() ? ObjSense::kMaximize : ObjSense::kMinimize;

    hm.lp_.col_cost_ = model.costs();
    hm.lp_.col_lower_.resize(n);
    hm.lp_.col_upper_.resize(n);
    for (std::size_t j = 0; j < n; ++j) {
        hm.lp_.col_lower_[j] = model.lower_is_free(j) ? -kHighsInf : model.lower(j);
        hm.lp_.col_upper_[j] = model.upper_is_free(j) ? kHighsInf : model.upper(j);
    }

    hm.lp_.row_lower_.resize(m);
    hm.lp_.row_upper_.resize(m);
    // HighsSparseMatrix carries its OWN dimensions and they are NOT inferred
    // from HighsLp. Leaving them at 0 makes HiGHS read an empty matrix, so
    // every equality row degenerates to 0 = rhs and the model comes back
    // Infeasible while the dense tableau solves it happily.
    hm.lp_.a_matrix_.num_col_ = static_cast<HighsInt>(n);
    hm.lp_.a_matrix_.num_row_ = static_cast<HighsInt>(m);
    hm.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    // HighsSparseMatrix DEFAULT-CONSTRUCTS with start_ = {0}. Appending the
    // leading zero without clearing yields start_ = {0, 0, nnz}, so HiGHS reads
    // row 0 as spanning [0,0) and reports "0 nonzeros" for the whole matrix --
    // every equality then degenerates to 0 = rhs and the model is Infeasible.
    // The symptom is a clean Infeasible on a model the dense tableau solves.
    hm.lp_.a_matrix_.start_.clear();
    hm.lp_.a_matrix_.start_.reserve(m + 1);
    hm.lp_.a_matrix_.index_.reserve(model.num_nonzeros());
    hm.lp_.a_matrix_.value_.reserve(model.num_nonzeros());
    hm.lp_.a_matrix_.start_.push_back(0);
    for (std::size_t i = 0; i < m; ++i) {
        const double b = model.rhs(i);
        switch (model.sense(i)) {
            case LpSense::LE:
                hm.lp_.row_lower_[i] = -kHighsInf;
                hm.lp_.row_upper_[i] = b;
                break;
            case LpSense::GE:
                hm.lp_.row_lower_[i] = b;
                hm.lp_.row_upper_[i] = kHighsInf;
                break;
            default:
                hm.lp_.row_lower_[i] = b;
                hm.lp_.row_upper_[i] = b;
                break;
        }
        for (std::size_t k = model.row_begin(i); k < model.row_end(i); ++k) {
            hm.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(model.col_at(k)));
            hm.lp_.a_matrix_.value_.push_back(model.val_at(k));
        }
        hm.lp_.a_matrix_.start_.push_back(static_cast<HighsInt>(hm.lp_.a_matrix_.index_.size()));
    }

    Highs highs;
    highs.setOptionValue("output_flag", false);
    // kWarning is BENIGN and common -- HiGHS raises it for things like an
    // unscaled model or a presolve remark, and the solve still returns an
    // optimal basis. Treating anything other than kOk as failure made every
    // feasible model here report Infeasible while the dense tableau solved it.
    // Only kError means the call did not happen.
    if (highs.passModel(hm) == HighsStatus::kError) {
        out.status = LpStatus::Infeasible;
        return out;
    }
    if (highs.run() == HighsStatus::kError) {
        out.status = LpStatus::Infeasible;
        return out;
    }

    const HighsModelStatus st = highs.getModelStatus();
    if (st == HighsModelStatus::kOptimal) {
        out.status = LpStatus::Optimal;
    } else if (st == HighsModelStatus::kUnbounded) {
        out.status = LpStatus::Unbounded;
        return out;
    } else if (st == HighsModelStatus::kIterationLimit ||
               st == HighsModelStatus::kTimeLimit) {
        out.status = LpStatus::IterationLimit;
        return out;
    } else {
        out.status = LpStatus::Infeasible;
        return out;
    }

    out.objective = highs.getInfo().objective_function_value;
    out.iterations = static_cast<std::size_t>(highs.getInfo().simplex_iteration_count);
    out.x = highs.getSolution().col_value;
    return out;
}

#endif  // LINE_MP_HAVE_HIGHS

/**
 * Solve, choosing the backend by arithmetic and size.
 *
 * Rational and every extended-precision T always take the exact dense path, at
 * compile time, because that is where this tree's exactness guarantees live.
 * A double model goes to HiGHS only once it is wider than `dense_max_cols`,
 * so small models keep bit-for-bit the answers their goldens were taken with
 * and the two backends stay comparable on exactly the instances the tests use.
 */
template <class T>
LpSolution<T> lp_solve(const LpModel<T>& model, std::size_t dense_max_cols = 512) {
    (void)dense_max_cols;
    return simplex_solve(model);
}

#ifdef LINE_MP_HAVE_HIGHS
template <>
inline LpSolution<double> lp_solve(const LpModel<double>& model, std::size_t dense_max_cols) {
    if (model.num_vars() > dense_max_cols) return highs_solve(model);
    return simplex_solve(model);
}
#endif

}  // namespace lp
}  // namespace line

#endif  // LINE_UTIL_LP_HIGHS_H
