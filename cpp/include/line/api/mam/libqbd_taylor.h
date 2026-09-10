/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The algorithm ported here is libQBD by S. Astaf'ev (IAMR Karelian Research
 * Centre RAS), BSD-3-Clause; the license ships in
 * matlab/lib/thirdparty/libQBD/LICENSE.
 */
#ifndef LINE_API_MAM_LIBQBD_TAYLOR_H
#define LINE_API_MAM_LIBQBD_TAYLOR_H

/**
 * Transient distribution of a level-independent-in-the-tail QBD by an adaptive
 * Taylor series (libQBD `QBD` + `TaylorSeriesAdaptive`).
 *
 * `solver_mam_ldqbd_transient.m` uses this for the INFINITE-buffer case, where
 * there is no finite generator to exponentiate. The method advances
 * pi(t) in steps of 1/|min diagonal|, and at each step sums the Taylor series
 *
 *     pi(t + h) = sum_k (h Q)^k / k! pi(t)
 *
 * in the UNIFORMIZED time h = 1/min_elem, truncating when the tail bound
 *
 *     ||d_k||_1 P(k+2, 2) e^2 2^-(k+1)
 *
 * falls below the requested error. `P(a,x)` is the regularized lower incomplete
 * gamma function, which is what makes this a genuine a-posteriori bound rather
 * than a heuristic cutoff, and it is the one special function the port had to
 * add (`gammainc_lower` below).
 *
 * THE STATE SPACE GROWS AS IT MUST. The distribution is a list of per-level row
 * vectors; multiplying by the generator can push mass one level higher, so
 * `mull_by_row_vector` appends a level whenever the new top carries any mass.
 * That is how an infinite buffer is handled without a truncation parameter:
 * the represented depth is whatever the elapsed time has actually reached.
 *
 * Levels above the last supplied one REPEAT it. `get_A_*` clamps the index, so
 * `add_final_level` defines the repeating block and every deeper level reuses
 * it, which is the QBD structure itself rather than an approximation.
 *
 * ARITHMETIC. Gated on transcendental: the truncation test evaluates an
 * incomplete gamma and an exponential.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Regularized lower incomplete gamma P(a, x), MATLAB's `gammainc(x, a, 'lower')`.
 *
 * Series below the transition point, continued fraction above it (Numerical
 * Recipes 6.2); evaluated in double because it is a TRUNCATION TEST, not a
 * returned quantity -- the distribution itself is accumulated in T.
 */
inline double gammainc_lower(double a, double x) {
    if (x < 0.0 || a <= 0.0) throw InputError("gammainc_lower: a must be positive and x >= 0");
    if (x == 0.0) return 0.0;
    const double gln = std::lgamma(a);
    if (x < a + 1.0) {
        double ap = a, del = 1.0 / a, sum = del;
        for (int n = 0; n < 1000; ++n) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
        }
        return sum * std::exp(-x + a * std::log(x) - gln);
    }
    const double tiny = 1e-300;
    double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
    for (int i = 1; i <= 1000; ++i) {
        const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::fabs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-16) break;
    }
    return 1.0 - std::exp(-x + a * std::log(x) - gln) * h;
}

/** A QBD's level blocks, as libQBD's `QBD` class holds them. */
template <class T>
class LibQbdProcess {
  public:
    /** Level zero from its upward block alone; the local block is the row-sum negative. */
    void add_zero_level(const Matrix<T>& Aplus) {
        if (!A0_.empty() || !Ap_.empty())
            throw InputError("LibQbdProcess: level zero already exists");
        Ap_.push_back(Aplus);
        A0_.push_back(diag_negrowsum(Aplus));
    }
    void add_zero_level(const Matrix<T>& A0, const Matrix<T>& Aplus) {
        if (!A0_.empty() || !Ap_.empty())
            throw InputError("LibQbdProcess: level zero already exists");
        A0_.push_back(A0);
        Ap_.push_back(Aplus);
    }
    /** A level from its down and up blocks; the local block closes the rows. */
    void add_level(const Matrix<T>& Aminus, const Matrix<T>& Aplus) {
        check_filled();
        Am_.push_back(Aminus);
        A0_.push_back(diag_negrowsum2(Aminus, Aplus));
        Ap_.push_back(Aplus);
    }
    void add_level(const Matrix<T>& Aminus, const Matrix<T>& A0, const Matrix<T>& Aplus) {
        check_filled();
        Am_.push_back(Aminus);
        A0_.push_back(A0);
        Ap_.push_back(Aplus);
    }
    /** The repeating level: its up block is the previous one, reused for ever. */
    void add_final_level(const Matrix<T>& Aminus) {
        check_filled();
        const Matrix<T> prevAp = Ap_.back();
        Am_.push_back(Aminus);
        A0_.push_back(diag_negrowsum2(Aminus, prevAp));
        Ap_.push_back(prevAp);
    }
    void add_final_level(const Matrix<T>& Aminus, const Matrix<T>& A0) {
        check_filled();
        const Matrix<T> prevAp = Ap_.back();
        Am_.push_back(Aminus);
        A0_.push_back(A0);
        Ap_.push_back(prevAp);
    }

    bool empty() const { return A0_.empty(); }

    /** Blocks at a level, with every level above the last one REPEATING it. */
    const Matrix<T>& A0(std::size_t level) const {
        return A0_[std::min(level, A0_.size() - 1)];
    }
    const Matrix<T>& Aplus(std::size_t level) const {
        return Ap_[std::min(level, Ap_.size() - 1)];
    }
    const Matrix<T>& Aminus(std::size_t level) const {
        if (level == 0) throw InputError("LibQbdProcess: A_minus at level zero is undefined");
        return Am_[std::min(level, Am_.size()) - 1];
    }
    /** The most negative diagonal entry over every local block. */
    T min_element() const {
        T v = num_traits<T>::from_int(0);
        for (const Matrix<T>& M : A0_)
            for (std::size_t i = 0; i < M.rows(); ++i)
                if (M(i, i) < v) v = M(i, i);
        return v;
    }

    /**
     * vec Q scaled by `cons`, where `vec` is one row vector per level.
     *
     * The result may be ONE LEVEL LONGER than the input: mass pushed above the
     * current top is kept whenever it is nonzero, which is what lets the
     * representation grow with the elapsed time instead of being truncated.
     */
    std::vector<std::vector<T>> mul_row(const std::vector<std::vector<T>>& vec,
                                        const T& cons) const {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t n = vec.size();
        if (n == 0) throw InputError("LibQbdProcess: an empty vector was passed");
        std::vector<std::vector<T>> res(n + 1);
        for (std::size_t j = 0; j <= n; ++j) {
            const std::size_t dim =
                (j < n) ? A0(j).cols() : Aplus(n - 1).cols();
            res[j].assign(dim, zero);
            if (j > 0) accumulate(res[j], vec[j - 1], Aplus(j - 1));
            if (j < n) accumulate(res[j], vec[j], A0(j));
            if (j + 1 < n) accumulate(res[j], vec[j + 1], Aminus(j + 1));
            for (T& v : res[j]) v *= cons;
        }
        // The reference keeps the new top level only when it carries mass, and
        // only in the n >= 3 branch; below that it always appends.
        if (n >= 3) {
            double nrm = 0.0;
            for (const T& v : res[n]) nrm += std::fabs(num_traits<T>::to_double(v));
            if (!(nrm > 0.0)) res.pop_back();
        }
        return res;
    }

  private:
    static Matrix<T> diag_negrowsum(const Matrix<T>& A) {
        const T zero = num_traits<T>::from_int(0);
        Matrix<T> D(A.rows(), A.rows(), zero);
        for (std::size_t i = 0; i < A.rows(); ++i) {
            T s = zero;
            for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j);
            D(i, i) = -s;
        }
        return D;
    }
    static Matrix<T> diag_negrowsum2(const Matrix<T>& A, const Matrix<T>& B) {
        const T zero = num_traits<T>::from_int(0);
        Matrix<T> D(A.rows(), A.rows(), zero);
        for (std::size_t i = 0; i < A.rows(); ++i) {
            T s = zero;
            for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j);
            for (std::size_t j = 0; j < B.cols(); ++j) s += B(i, j);
            D(i, i) = -s;
        }
        return D;
    }
    static void accumulate(std::vector<T>& out, const std::vector<T>& v, const Matrix<T>& M) {
        if (v.size() != M.rows()) return;  // a level whose width does not meet this block
        for (std::size_t i = 0; i < M.rows(); ++i) {
            if (v[i] == num_traits<T>::from_int(0)) continue;
            for (std::size_t j = 0; j < M.cols() && j < out.size(); ++j) out[j] += v[i] * M(i, j);
        }
    }

    void check_filled() const {
        if (A0_.size() != Ap_.size() || Ap_.size() != Am_.size() + 1)
            throw InputError("LibQbdProcess: unfilled levels found");
    }

    std::vector<Matrix<T>> Ap_, A0_, Am_;
};

/** What the adaptive Taylor series returns: the reference grid and its laws. */
template <class T>
struct TaylorSeriesResult {
    std::vector<double> times;                            ///< the reference points
    std::vector<std::vector<std::vector<T>>> dists;       ///< per point, per level, per phase
};

/**
 * libQBD's `TaylorSeriesAdaptive`, restricted to the reference grid that
 * `solver_mam_ldqbd_transient` reads (`get_reference_times` and
 * `get_reference_dists`); the interpolation to arbitrary points is not ported
 * because no caller in this tree asks for it.
 *
 * @param pi0       initial law, one row vector per level
 * @param error     per-step truncation target (the reference passes options.tol)
 * @param max_time  advance until the grid covers this horizon
 * @param proc the level-dependent QBD being integrated
 */
template <class T>
TaylorSeriesResult<T> taylor_series_adaptive(const LibQbdProcess<T>& proc,
                                             const std::vector<std::vector<T>>& pi0,
                                             double error, double max_time) {
    static_assert(num_traits<T>::has_transcendental,
                  "taylor_series_adaptive evaluates an incomplete gamma truncation bound");
    if (proc.empty()) throw InputError("taylor_series_adaptive: the generator is empty");
    const double min_elem = -num_traits<T>::to_double(proc.min_element());
    if (!(min_elem > 0.0))
        throw NumericError("taylor_series_adaptive: the generator has no negative diagonal");
    const unsigned max_degree = 177u;  // libQBD's get_max_factor() for double

    TaylorSeriesResult<T> out;
    out.times.push_back(0.0);
    out.dists.push_back(pi0);

    const T min_elem_inv = num_traits<T>::from_double(1.0 / min_elem);
    while (out.times.back() < max_time) {
        std::vector<std::vector<T>> deriv = out.dists.back();
        std::vector<std::vector<T>> res = deriv;
        unsigned k = 0;
        double two_delta_in_n = 0.5;
        double er = std::numeric_limits<double>::infinity();
        while (er > error && k < max_degree) {
            deriv = proc.mul_row(deriv, min_elem_inv);
            // res += deriv / (k+1)!
            const double c = std::exp(-std::lgamma(static_cast<double>(k) + 2.0));
            const T ct = num_traits<T>::from_double(c);
            if (res.size() < deriv.size()) res.resize(deriv.size());
            for (std::size_t l = 0; l < deriv.size(); ++l) {
                if (res[l].size() < deriv[l].size())
                    res[l].resize(deriv[l].size(), num_traits<T>::from_int(0));
                for (std::size_t i = 0; i < deriv[l].size(); ++i) res[l][i] += ct * deriv[l][i];
            }
            double nrm = 0.0;
            for (const std::vector<T>& lv : deriv)
                for (const T& v : lv) nrm += std::fabs(num_traits<T>::to_double(v));
            er = nrm * gammainc_lower(static_cast<double>(k) + 2.0, 2.0) * std::exp(2.0) *
                 two_delta_in_n;
            two_delta_in_n *= 0.5;
            ++k;
        }
        out.dists.push_back(res);
        out.times.push_back(out.times.back() + 1.0 / min_elem);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_LIBQBD_TAYLOR_H
