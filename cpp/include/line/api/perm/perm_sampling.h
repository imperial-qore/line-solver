/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PERM_PERM_SAMPLING_H
#define LINE_API_PERM_PERM_SAMPLING_H

/**
 * RANDOMIZED permanents: the AdaPart rejection sampler and the Huber-Law
 * acceptance-rejection importance sampler.
 *
 * Port of python/line_solver/api/perm/sampling.py, itself a twin of
 * jline.lib.perm.AdaPartSampler and jline.lib.perm.HuberLawSampler, and of
 * MATLAB's perm_adapart.m and perm_huberlaw.m.
 *
 * WHAT THESE ARE FOR. The exact algorithms in permanent.h cost 2^n or n!, and
 * the two schemes in perm_approx.h buy their speed with a bias no caller can
 * bound. These two are unbiased Monte Carlo estimators instead: the estimate
 * carries sampling noise, but no systematic error, so averaging more draws
 * converges to the permanent rather than to a nearby number.
 *
 *  - `perm_adapart` recursively partitions the permutation space, bounds each
 *    part by the Soules column bound, and draws a part with probability
 *    proportional to its bound. The acceptance ratio scales the root bound into
 *    an unbiased estimate.
 *  - `perm_huberlaw` rescales the matrix to doubly stochastic, draws a
 *    permutation column by column under the Huber-Law bound on the remaining
 *    permanent, and multiplies the acceptance ratio by the rescaling constant.
 *
 * TWO ADAPART FIXES ARE CARRIED HERE, and both are load-bearing rather than
 * cosmetic:
 *
 *  1. The column search considers ONLY the columns still unassigned. Scoring an
 *     already-assigned column just re-derives its own constraint, which always
 *     looks cheapest, so the sampler re-splits the same column forever and never
 *     completes an assignment.
 *  2. The expansion discounts the bound of the element ACTUALLY REMOVED, not the
 *     root bound. With the root bound the running total drifts and the "refine
 *     until improved" loop cannot exit.
 *
 * Even so the Soules bound is TIGHT on a matrix with equal entries, so no
 * refinement can improve it; the loop therefore stops on the first
 * non-improving expansion, which is exactly the shape a replicated demand
 * matrix produces (see pfqn_jointmarg).
 *
 * BOTH REQUIRE A NONNEGATIVE MATRIX, and both need FULL SUPPORT: an exact zero
 * makes the Huber-Law rescaling degenerate and starves the AdaPart acceptance.
 * The callers in the pfqn layer refuse a structurally zero matrix rather than
 * flooring it.
 *
 * REPRODUCIBILITY: the draws come from a seeded std::mt19937_64, so a run with
 * a fixed seed repeats. The Python twin draws from a numpy Generator and the
 * JAR from java.util.Random, so the three agree in distribution but not sample
 * by sample.
 *
 * ARITHMETIC: double. Both are floating-point Monte Carlo schemes.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace perm {

namespace samplingdetail {

/** n! in double by the plain product, as in the JAR and the Python twin. */
inline double factorial_plain(std::size_t n) {
    double result = 1.0;
    for (std::size_t i = 1; i <= n; ++i) result *= static_cast<double>(i);
    return result;
}

/** Rejects a matrix that is not square or carries a negative entry. */
inline void require_square_nonnegative(const Matrix<double>& m, const char* who) {
    if (m.rows() != m.cols())
        throw InputError(std::string(who) + ": the matrix must be square, got " +
                         std::to_string(m.rows()) + " by " + std::to_string(m.cols()));
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            if (m(i, j) < 0.0)
                throw InputError(std::string(who) + ": entry (" + std::to_string(i + 1) + "," +
                                 std::to_string(j + 1) + ") is negative");
}

/**
 * Refuse a matrix the samplers cannot take. Twin of the approx-header guard.
 *
 * A zero used to be floored before scaling, and that substitution is not
 * invertible: it changes the permanent by n!*eps, which is O(1) by n = 18.
 * Positivity is sufficient but not necessary -- the sharp precondition is
 * TOTAL SUPPORT -- but it is O(n^2) and is the contract the headers state.
 */
inline void require_full_support(const Matrix<double>& m, const char* who) {
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            if (!(m(i, j) > 0.0))
                throw InputError(std::string(who) +
                                 ": requires a strictly positive matrix, but entry (" +
                                 std::to_string(i + 1) + "," + std::to_string(j + 1) +
                                 ") is " + std::to_string(m(i, j)) +
                                 ", so the matrix has no full support. Flooring it would change"
                                 " the permanent by n!*eps, which is O(1) by n=18. Use the exact"
                                 " engine.");
}

/**
 * Maximum-weight perfect assignment, by the O(n^3) Hungarian algorithm.
 *
 * Replaces the row-by-row greedy that used to stand in for it. The greedy
 * returns a zero-weight assignment on inputs that admit a positive one: on
 * [[1,2],[0,3]] row 0 takes the larger entry in column 1, leaving row 1 with
 * the zero. That weight is alpha3, which sets the flooring level
 * alpha1 = alpha3*delta/(3 n!) of the Huber-Law bound, so a suboptimal
 * assignment weakens the method's own guarantee. The MATLAB, JAR and python
 * twins solve the same problem and agree on the optimal VALUE, which is all
 * alpha3 depends on; they need not agree on the permutation when it is
 * degenerate.
 */
inline std::vector<std::size_t> max_weight_assignment(const Matrix<double>& weight) {
    const std::size_t n = weight.rows();
    std::vector<std::size_t> assignment(n, 0);
    if (n == 0) return assignment;

    const double kInf = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> cost(n + 1, std::vector<double>(n + 1, 0.0));
    for (std::size_t i = 1; i <= n; ++i)
        for (std::size_t j = 1; j <= n; ++j) cost[i][j] = -weight(i - 1, j - 1);

    std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0);
    std::vector<std::size_t> p(n + 1, 0), way(n + 1, 0);

    for (std::size_t i = 1; i <= n; ++i) {
        p[0] = i;
        std::size_t j0 = 0;
        std::vector<double> minv(n + 1, kInf);
        std::vector<bool> used(n + 1, false);
        do {
            used[j0] = true;
            const std::size_t i0 = p[j0];
            std::size_t j1 = 0;
            double delta = kInf;
            for (std::size_t j = 1; j <= n; ++j) {
                if (!used[j]) {
                    const double cur = cost[i0][j] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for (std::size_t j = 0; j <= n; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            const std::size_t j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    for (std::size_t j = 1; j <= n; ++j)
        if (p[j] != 0) assignment[p[j] - 1] = j - 1;
    return assignment;
}

/**
 * Soules upper bound of the permanent, a product of column bounds.
 *
 * gamma(k) = (k!)^(1/k), delta(i) = gamma(n-i) - gamma(n-i-1), and each column
 * is sorted ascending before the weights are applied.
 */
inline double soules_bound(const Matrix<double>& m) {
    const std::size_t n = m.rows();
    if (n == 0) return 1.0;
    std::vector<double> gamma(n + 1, 0.0);
    double fact = 1.0;
    for (std::size_t k = 1; k <= n; ++k) {
        fact *= static_cast<double>(k);
        gamma[k] = std::pow(fact, 1.0 / static_cast<double>(k));
    }
    std::vector<double> delta(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) delta[i] = gamma[n - i] - gamma[n - i - 1];

    double prod = 1.0;
    std::vector<double> col(n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) col[i] = m(i, j);
        std::sort(col.begin(), col.end());
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) s += delta[i] * col[i];
        prod *= s;
    }
    return prod;
}

}  // namespace samplingdetail

/**
 * Adaptive partitioning (AdaPart) sampler for the permanent.
 *
 * The partial assignment of the partition is a vector t of length n whose entry
 * j is the row assigned to column j, or n when column j is still free.
 */
class AdaPartSampler {
  public:
    /** Draw budgets; 'classic' is the reference default. */
    enum class Mode { Classic, Time, Sample };

    /**
     * @param matrix                    nonnegative square matrix
     * @param maximum_accepted_samples  acceptance budget of Classic
     * @param maximum_time              time budget in milliseconds of Time
     * @param maximum_samples           draw budget of Sample
     * @param mode                      which budget applies
     * @param seed                      seed of the draws
     */
    explicit AdaPartSampler(const Matrix<double>& matrix, int maximum_accepted_samples = 100,
                            double maximum_time = 30000.0, int maximum_samples = 450,
                            Mode mode = Mode::Classic, std::uint64_t seed = 0)
        : matrix_(matrix),
          n_(matrix.rows()),
          maximum_accepted_samples_(maximum_accepted_samples),
          maximum_time_(maximum_time),
          maximum_samples_(maximum_samples),
          mode_(mode),
          rng_(seed),
          value_(0.0) {
        samplingdetail::require_square_nonnegative(matrix_, "AdaPartSampler");
        if (n_ > 0) samplingdetail::require_full_support(matrix_, "AdaPartSampler");
    }

    /** Run the sampler in the configured mode and return the estimate. */
    double solve() {
        const double z_ub = samplingdetail::soules_bound(matrix_);
        long accepted = 0, total = 0;
        const std::clock_t start = std::clock();
        // Bounded independently of the scaling: with perm(A) = 0 the acceptance
        // probability is 0 and Classic mode would never terminate. A cap that
        // RETURNS a number would be a workaround, so it throws.
        const long kMaxDraws = 1000000;
        while (true) {
            if (mode_ == Mode::Classic && accepted >= maximum_accepted_samples_) break;
            if (mode_ == Mode::Classic && total >= kMaxDraws)
                throw InputError("perm_adapart: only " + std::to_string(accepted) + " of the " +
                                 std::to_string(maximum_accepted_samples_) +
                                 " required acceptances were obtained in " +
                                 std::to_string(total) +
                                 " draws. Use the exact engine.");
            if (mode_ == Mode::Sample && total >= maximum_samples_) break;
            if (mode_ == Mode::Time && elapsed_ms(start) >= maximum_time_) break;
            accepted += sample(start);
            ++total;
        }
        value_ = (total > 0) ? z_ub * static_cast<double>(accepted) / static_cast<double>(total)
                             : 0.0;
        return value_;
    }

    /** Estimate of the last solve. */
    double value() const { return value_; }

  private:
    typedef std::vector<std::size_t> Assignment;

    static double elapsed_ms(std::clock_t start) {
        return 1000.0 * static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;
    }

    bool within_time(std::clock_t start) const {
        if (mode_ != Mode::Time) return true;
        return elapsed_ms(start) < maximum_time_;
    }

    /** True while some partition element still has an unassigned column. */
    bool any_free(const std::set<Assignment>& s_set) const {
        for (std::set<Assignment>::const_iterator it = s_set.begin(); it != s_set.end(); ++it)
            if (has_free(*it)) return true;
        return false;
    }

    bool has_free(const Assignment& t) const {
        for (std::size_t j = 0; j < t.size(); ++j)
            if (t[j] == n_) return true;
        return false;
    }

    /** Zero out the entries excluded by the partial assignment t. */
    Matrix<double> modify_matrix(const Matrix<double>& m, const Assignment& t) const {
        Matrix<double> out(n_, n_, 0.0);
        std::vector<bool> row_used(n_, false);
        for (std::size_t j = 0; j < t.size(); ++j)
            if (t[j] != n_) row_used[t[j]] = true;
        for (std::size_t j = 0; j < t.size(); ++j) {
            if (t[j] != n_) {
                out(t[j], j) = m(t[j], j);
            } else {
                for (std::size_t i = 0; i < n_; ++i)
                    if (!row_used[i]) out(i, j) = m(i, j);
            }
        }
        return out;
    }

    /**
     * Pick the column whose expansion minimizes the summed Soules bound.
     *
     * Only columns still unassigned in s_sub are candidates; see the header.
     */
    std::size_t select_column(const Matrix<double>& s_matrix, double removed_ub, double ub,
                              const Assignment& s_sub, double* new_ub) const {
        double best = std::numeric_limits<double>::infinity();
        std::size_t best_col = 0;
        bool found = false;
        for (std::size_t i = 0; i < n_; ++i) {
            if (s_sub[i] != n_) continue;
            double total = 0.0;
            for (std::size_t j = 0; j < n_; ++j) {
                Assignment a(n_, n_);
                a[i] = j;
                total += samplingdetail::soules_bound(modify_matrix(s_matrix, a));
            }
            if (!found || total < best) {
                best = total;
                best_col = i;
                found = true;
            }
        }
        if (!found) {
            *new_ub = ub;
            return 0;
        }
        *new_ub = ub - removed_ub + best;
        return best_col;
    }

    /** Draw a partition element, or the slack index that means rejection. */
    std::size_t compute_probabilities(const std::set<Assignment>& s_set, double zub_s) {
        std::vector<double> p;
        p.reserve(s_set.size() + 1);
        double sum = 0.0;
        for (std::set<Assignment>::const_iterator it = s_set.begin(); it != s_set.end(); ++it) {
            const double b = samplingdetail::soules_bound(modify_matrix(matrix_, *it));
            p.push_back(b);
            sum += b;
        }
        if (sum > 0.0 && zub_s > 0.0) {
            double norm = 0.0;
            for (std::size_t i = 0; i < p.size(); ++i) {
                p[i] /= zub_s;
                norm += p[i];
            }
            p.push_back(1.0 - norm);
        } else {
            p.push_back(1.0);
        }
        double total = 0.0;
        for (std::size_t i = 0; i < p.size(); ++i) total += std::fabs(p[i]);
        if (total > 0.0)
            for (std::size_t i = 0; i < p.size(); ++i) p[i] = std::fabs(p[i]) / total;

        const double u = uniform_();
        double cum = 0.0;
        for (std::size_t i = 0; i < p.size(); ++i) {
            cum += p[i];
            if (u <= cum) return i;
        }
        return p.size() - 1;
    }

    /** Keep the drawn element, completing it when one column is left. */
    std::set<Assignment> subset(const std::set<Assignment>& s_set, std::size_t c) const {
        std::set<Assignment>::const_iterator it = s_set.begin();
        std::advance(it, static_cast<long>(c));
        Assignment s_inter = *it;
        std::size_t free_count = 0, free_pos = 0;
        for (std::size_t j = 0; j < s_inter.size(); ++j)
            if (s_inter[j] == n_) {
                ++free_count;
                free_pos = j;
            }
        if (free_count == 1) {
            std::vector<bool> used(n_, false);
            for (std::size_t j = 0; j < s_inter.size(); ++j)
                if (s_inter[j] != n_) used[s_inter[j]] = true;
            for (std::size_t v = 0; v < n_; ++v)
                if (!used[v]) {
                    s_inter[free_pos] = v;
                    break;
                }
        }
        std::set<Assignment> out;
        out.insert(s_inter);
        return out;
    }

    /** Draw one partition path, returning 1 if accepted and 0 if rejected. */
    int sample(std::clock_t start) {
        std::set<Assignment> s_set;
        s_set.insert(Assignment(n_, n_));
        while (any_free(s_set) && within_time(start)) {
            const Assignment s_init = *s_set.begin();
            const double zub_s = samplingdetail::soules_bound(modify_matrix(matrix_, s_init));
            double ub = zub_s;
            bool init = true;
            while ((ub >= zub_s || init) && within_time(start)) {
                init = false;
                // Only elements with a free position can be refined; expanding a
                // complete assignment yields no children and stalls the sampler.
                std::vector<Assignment> s_list;
                for (std::set<Assignment>::const_iterator it = s_set.begin(); it != s_set.end();
                     ++it)
                    if (has_free(*it)) s_list.push_back(*it);
                if (s_list.empty()) break;
                const std::size_t pick =
                    static_cast<std::size_t>(uniform_() * static_cast<double>(s_list.size()));
                const Assignment s_sub = s_list[pick < s_list.size() ? pick : s_list.size() - 1];
                s_set.erase(s_sub);
                const Matrix<double> sub_matrix = modify_matrix(matrix_, s_sub);
                // Discount the bound of the element actually removed, not the
                // root bound; see the header.
                const double sub_ub = samplingdetail::soules_bound(sub_matrix);
                double new_ub = ub;
                const std::size_t j = select_column(sub_matrix, sub_ub, ub, s_sub, &new_ub);
                for (std::size_t i = 0; i < n_; ++i) {
                    bool taken = false;
                    for (std::size_t k = 0; k < s_sub.size(); ++k)
                        if (s_sub[k] == i) taken = true;
                    if (taken) continue;
                    Assignment s_add = s_sub;
                    s_add[j] = i;
                    s_set.insert(s_add);
                }
                const bool no_progress = new_ub >= ub;
                ub = new_ub;
                // The Soules bound is tight on matrices with equal entries, so
                // refinement cannot improve it and "refine until improved" would
                // never exit. Stop on the first non-improving expansion instead;
                // a tight bound means the draw is accepted with probability 1.
                if (no_progress) break;
            }
            const std::size_t c = compute_probabilities(s_set, zub_s);
            if (c == s_set.size()) return 0;
            s_set = subset(s_set, c);
        }
        return 1;
    }

    double uniform_() {
        return std::generate_canonical<double, 53>(rng_);
    }

    Matrix<double> matrix_;
    std::size_t n_;
    int maximum_accepted_samples_;
    double maximum_time_;
    int maximum_samples_;
    Mode mode_;
    std::mt19937_64 rng_;
    double value_;
};

/**
 * Huber-Law acceptance-rejection sampler for the permanent.
 *
 * The matrix is rescaled to be doubly stochastic, a permutation is drawn column
 * by column under the Huber-Law upper bound on the remaining permanent, and the
 * acceptance ratio times the rescaling constant estimates the permanent.
 */
class HuberLawSampler {
  public:
    enum class Mode { Classic, Time, Sample };

    /**
     * @param matrix             nonnegative square matrix
     * @param delta              relative accuracy target, sets the budget K
     * @param alpha2             convergence threshold of the rescaling
     * @param epsilon            failure probability target, sets the budget K
     * @param mode               which budget applies
     * @param number_of_samples  draw budget of Sample
     * @param maximum_time       time budget in milliseconds of Time
     * @param seed               seed of the draws
     */
    explicit HuberLawSampler(const Matrix<double>& matrix, double delta = 0.1,
                             double alpha2 = 0.000001, double epsilon = 0.1,
                             Mode mode = Mode::Classic, int number_of_samples = 1000,
                             double maximum_time = 30000.0, std::uint64_t seed = 0)
        : matrix_(matrix),
          n_(matrix.rows()),
          delta_(delta),
          alpha2_(alpha2),
          epsilon_(epsilon),
          mode_(mode),
          number_of_samples_(number_of_samples),
          maximum_time_(maximum_time),
          rng_(seed),
          c_matrix_(matrix.rows(), matrix.cols(), 0.0),
          rescaling_constant_(1.0),
          value_(0.0) {
        samplingdetail::require_square_nonnegative(matrix_, "HuberLawSampler");
        if (n_ > 0) samplingdetail::require_full_support(matrix_, "HuberLawSampler");
    }

    /** Run the sampler in the configured mode and return the estimate. */
    double solve() {
        rescale();
        const long k = static_cast<long>(14.0 * std::pow(delta_, -2.0) *
                                         std::log(2.0 / epsilon_));
        const std::clock_t start = std::clock();
        long accepted = 0, total = 0;
        // Bounded independently of the scaling, as in the AdaPart sampler.
        const long kMaxDraws = 1000000;
        while (true) {
            if (mode_ == Mode::Classic && accepted >= k) break;
            if (mode_ == Mode::Classic && total >= kMaxDraws)
                throw InputError("perm_huberlaw: only " + std::to_string(accepted) + " of the " +
                                 std::to_string(k) +
                                 " required acceptances were obtained in " +
                                 std::to_string(total) +
                                 " draws. Relax delta or use the exact engine.");
            if (mode_ == Mode::Sample && total >= number_of_samples_) break;
            if (mode_ == Mode::Time && elapsed_ms(start) >= maximum_time_) break;
            accepted += sample();
            ++total;
        }
        value_ = (total > 0) ? static_cast<double>(accepted) / static_cast<double>(total) *
                                   rescaling_constant_
                             : 0.0;
        return value_;
    }

    /** Estimate of the last solve. */
    double value() const { return value_; }

  private:
    static double elapsed_ms(std::clock_t start) {
        return 1000.0 * static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;
    }

    /** Huber-Law bound factor of a row with remaining mass r. */
    static double h(double r) {
        if (r >= 1.0) return r + 0.5 * std::log(r) + M_E - 1.0;
        return 1.0 + (M_E - 1.0) * r;
    }

    /** Unnormalized selection weights of each row for column j. */
    std::vector<double> precomputing(const Matrix<double>& m, std::size_t j) const {
        std::vector<double> hr(n_, 0.0), c(n_, 0.0);
        double hr_product = 1.0;
        for (std::size_t i = 0; i < n_; ++i) {
            c[i] = m(i, j);
            double rowsum = 0.0;
            for (std::size_t k = 0; k < n_; ++k) rowsum += m(i, k);
            hr[i] = h(rowsum - c[i]);
            hr_product *= hr[i];
        }
        const double exp_factor = std::exp(static_cast<double>(n_) - 1.0);
        std::vector<double> out(n_, 0.0);
        for (std::size_t i = 0; i < n_; ++i)
            out[i] = (hr[i] != 0.0) ? hr_product / hr[i] * c[i] / exp_factor : 0.0;
        return out;
    }

    /** Draw one permutation; a rejected draw is reported by the return value. */
    int sample() {
        Matrix<double> m = c_matrix_;
        for (std::size_t j = 0; j < n_; ++j) {
            std::vector<double> p = precomputing(m, j);
            double rowprod = 1.0;
            for (std::size_t i = 0; i < n_; ++i) {
                double rowsum = 0.0;
                for (std::size_t k = 0; k < n_; ++k) rowsum += m(i, k);
                rowprod *= h(rowsum);
            }
            const double ub = rowprod / std::exp(static_cast<double>(n_));
            std::vector<double> prob(n_ + 1, 0.0);
            double psum = 0.0;
            for (std::size_t i = 0; i < n_; ++i) {
                prob[i] = (ub > 0.0) ? p[i] / ub : 0.0;
                psum += prob[i];
            }
            prob[n_] = 1.0 - psum;
            if (prob[n_] < 0.0) {
                if (psum > 0.0)
                    for (std::size_t i = 0; i < n_; ++i) prob[i] /= psum;
                prob[n_] = 0.0;
            }
            const double u = uniform_();
            std::size_t selected = n_;
            double cum = 0.0;
            for (std::size_t i = 0; i <= n_; ++i) {
                cum += prob[i];
                if (u <= cum) {
                    selected = i;
                    break;
                }
            }
            if (selected == n_) return 0;  // rejected
            Matrix<double> next(n_, n_, 0.0);
            for (std::size_t a = 0; a < n_; ++a)
                for (std::size_t b = 0; b < n_; ++b)
                    if (a != selected && b != j) next(a, b) = m(a, b);
            next(selected, j) = m(selected, j);
            m = next;
        }
        return 1;
    }

    /** Greedy row-by-row assignment maximizing the cost, as in the JAR. */
    /** Alternate column and row normalization, filling the scaling diagonals. */
    Matrix<double> make_doubly_stochastic(const Matrix<double>& m, std::vector<double>* x,
                                          std::vector<double>* y) const {
        Matrix<double> result = m;
        x->assign(n_, 1.0);
        y->assign(n_, 1.0);
        double max_row_error = std::numeric_limits<double>::infinity();
        double max_col_error = std::numeric_limits<double>::infinity();
        // Capped: a row that sums to zero leaves max_row_error at 1 forever and
        // the guarded normalization below skips it, so this loop used to spin
        // without terminating. A cap that RETURNS is a workaround; this throws.
        const std::size_t kMaxSweeps = 10000;
        std::size_t sweeps = 0;
        while (max_row_error > alpha2_ || max_col_error > alpha2_) {
            if (++sweeps > kMaxSweeps)
                throw InputError(
                    "make_doubly_stochastic: did not converge in " +
                    std::to_string(kMaxSweeps) + " sweeps (row error " +
                    std::to_string(max_row_error) + ", column error " +
                    std::to_string(max_col_error) + " against a tolerance of " +
                    std::to_string(alpha2_) +
                    "). The usual cause is a matrix without total support.");
            for (std::size_t j = 0; j < n_; ++j) {
                double s = 0.0;
                for (std::size_t i = 0; i < n_; ++i) s += result(i, j);
                if (s > 0.0) {
                    for (std::size_t i = 0; i < n_; ++i) result(i, j) /= s;
                    (*y)[j] /= s;
                }
            }
            for (std::size_t i = 0; i < n_; ++i) {
                double s = 0.0;
                for (std::size_t j = 0; j < n_; ++j) s += result(i, j);
                if (s > 0.0) {
                    for (std::size_t j = 0; j < n_; ++j) result(i, j) /= s;
                    (*x)[i] /= s;
                }
            }
            max_col_error = 0.0;
            max_row_error = 0.0;
            for (std::size_t j = 0; j < n_; ++j) {
                double s = 0.0;
                for (std::size_t i = 0; i < n_; ++i) s += result(i, j);
                max_col_error = std::max(max_col_error, std::fabs(s - 1.0));
            }
            for (std::size_t i = 0; i < n_; ++i) {
                double s = 0.0;
                for (std::size_t j = 0; j < n_; ++j) s += result(i, j);
                max_row_error = std::max(max_row_error, std::fabs(s - 1.0));
            }
        }
        return result;
    }

    /** Rescale to doubly stochastic and set the scaling constant. */
    void rescale() {
        Matrix<double> log_matrix(n_, n_, 0.0);
        double max_element = 0.0;
        for (std::size_t i = 0; i < n_; ++i)
            for (std::size_t j = 0; j < n_; ++j) {
                // strictly positive here (require_full_support), so no floor
                log_matrix(i, j) = std::log(matrix_(i, j));
                max_element = std::max(max_element, matrix_(i, j));
            }
        if (!(max_element > 0.0))
            throw InputError("HuberLawSampler: the matrix is identically zero");

        const std::vector<std::size_t> assignment =
            samplingdetail::max_weight_assignment(log_matrix);

        Matrix<double> m_scaled(n_, n_, 0.0);
        for (std::size_t i = 0; i < n_; ++i)
            for (std::size_t j = 0; j < n_; ++j) m_scaled(i, j) = matrix_(i, j) / max_element;

        // alpha3 is a permanent lower bound of the SCALED matrix, the one floored below
        double alpha3 = 1.0;
        for (std::size_t i = 0; i < n_; ++i) alpha3 *= m_scaled(i, assignment[i]);
        const double alpha1 = alpha3 * delta_ / 3.0 / samplingdetail::factorial_plain(n_);
        for (std::size_t i = 0; i < n_; ++i)
            for (std::size_t j = 0; j < n_; ++j)
                m_scaled(i, j) = std::max(m_scaled(i, j), alpha1);

        std::vector<double> x, y;
        const Matrix<double> ds = make_doubly_stochastic(m_scaled, &x, &y);

        std::vector<double> z(n_, 1.0);
        for (std::size_t i = 0; i < n_; ++i) {
            double rowmax = 0.0;
            for (std::size_t j = 0; j < n_; ++j) rowmax = std::max(rowmax, ds(i, j));
            z[i] = (rowmax > 0.0) ? 1.0 / rowmax : 1.0;
        }
        for (std::size_t i = 0; i < n_; ++i)
            for (std::size_t j = 0; j < n_; ++j) c_matrix_(i, j) = z[i] * ds(i, j);

        double h_product = 1.0;
        for (std::size_t i = 0; i < n_; ++i) {
            double rowsum = 0.0;
            for (std::size_t j = 0; j < n_; ++j) rowsum += c_matrix_(i, j);
            h_product *= h(rowsum) / M_E;
        }
        double diagonal_product = 1.0;
        for (std::size_t i = 0; i < n_; ++i) diagonal_product *= x[i] * y[i] * z[i];
        rescaling_constant_ = h_product / diagonal_product *
                              std::pow(max_element, static_cast<double>(n_));
    }

    double uniform_() {
        return std::generate_canonical<double, 53>(rng_);
    }

    Matrix<double> matrix_;
    std::size_t n_;
    double delta_;
    double alpha2_;
    double epsilon_;
    Mode mode_;
    int number_of_samples_;
    double maximum_time_;
    std::mt19937_64 rng_;
    Matrix<double> c_matrix_;
    double rescaling_constant_;
    double value_;
};

/**
 * AdaPart estimate of the permanent, with the reference defaults.
 *
 * Twin of MATLAB perm_adapart.m and of python perm/sampling.py.
 */
inline double perm_adapart(const Matrix<double>& m, std::uint64_t seed = 0) {
    if (m.rows() == 0) return 1.0;
    AdaPartSampler s(m, 100, 30000.0, 450, AdaPartSampler::Mode::Classic, seed);
    return s.solve();
}

/**
 * Huber-Law estimate of the permanent, with the reference defaults.
 *
 * Twin of MATLAB perm_huberlaw.m and of python perm/sampling.py.
 */
inline double perm_huberlaw(const Matrix<double>& m, std::uint64_t seed = 0) {
    if (m.rows() == 0) return 1.0;
    HuberLawSampler s(m, 0.1, 0.000001, 0.1, HuberLawSampler::Mode::Classic, 1000, 30000.0, seed);
    return s.solve();
}

}  // namespace perm
}  // namespace line

#endif  // LINE_API_PERM_PERM_SAMPLING_H
