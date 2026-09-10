/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_TTL_LRUA_H
#define LINE_API_CACHE_TTL_LRUA_H

/**
 * TTL (characteristic-time) approximation of an LRU cache whose lists form an
 * arbitrary access graph.
 *
 * Templated port of matlab/src/api/cache/cache_ttl_lrua.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_ttl_lrua.java.
 *
 * Each item moves over the h+1 nodes "not cached" (node 0) and "in list l"
 * (node l), driven by the access graph R and by exponential timer races: given
 * characteristic times x_1..x_h, an item in list l is promoted along R with
 * probability 1 - exp(-lambda_l x_l) and demoted with probability
 * exp(-lambda_l x_l). The embedded chain over the reachable nodes is solved by
 * dtmc_solve, the mean holding times are 1/lambda_0 at node 0 and
 * (1 - exp(-lambda_l x_l))/lambda_l in list l, and the time-stationary
 * probabilities are the holding-time weighted normalization of the two. The
 * times are then fixed by sum_k prob(k,l) = m_l, one equation per list.
 *
 * HOW THE SYSTEM IS SOLVED: as in cache_t_lrum_map.h and cache_t_hlru.h, the
 * occupancy of list l is increasing in its own characteristic time, so the
 * system is h scalar bracketed root problems swept Gauss-Seidel, closed by
 * bisection. This is deterministic. MATLAB instead calls fsolve from a RANDOM
 * initial point -- rng(seed,'twister') with a default seed of 23000 and
 * x = 10*rand(1,h) -- so its answer depends on the seed argument and on the
 * Optimization Toolbox, and it stops at fsolve's default tolerance (the
 * reference instance in the tests leaves capacity residuals of ~1.4e-7). The
 * JAR uses a damped Newton on log(T) with a finite-difference Jacobian, which
 * is deterministic but still needs a Jacobian and a line search.
 *
 * ARITHMETIC: exp and a bisection tolerance, so transcendental arithmetic is
 * required.
 *
 * LATENT CONSTRAINT IN THE REFERENCE: MATLAB reads the demotion probability as
 * exp(-lambda(1,i,k)*x(k-1)) for every k that is a successor of some node j,
 * which for k = 1 (the "not cached" node, 0 here) indexes x(0) and is a hard
 * MATLAB error. An access graph that routes back into the not-cached node is
 * therefore not expressible; this port raises InputError on it instead of
 * relying on an index-out-of-range.
 *
 * ONE-USER REDUCTION: the MATLAB signature takes a (u x n x h+1) array but its
 * body reads lambda(1,i,j) only, so every user beyond the first is ignored.
 * The port takes the (n x h+1) matrix that MATLAB actually uses, which makes
 * the reduction explicit rather than silent. The caller in
 * solver_mva_cache_analyzer.m passes per-user rates whose relevant slice is
 * already the aggregate.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/rootfind.h"

namespace line {
namespace cache {

namespace detail {

/**
 * Time-stationary probability of each item at each node, MATLAB's randprob,
 * for the characteristic times x.
 */
template <class T>
Matrix<T> lrua_randprob(const Matrix<T>& lambda, const std::vector<Matrix<T>>& R,
                        const std::vector<T>& x) {
    using std::exp;
    const std::size_t n = lambda.rows();
    const std::size_t h = x.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    Matrix<T> randprob(n, h + 1, zero);

    for (std::size_t i = 0; i < n; ++i) {
        // no-arrivals item rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
        if (lambda(i, 0) == zero) {
            randprob(i, 0) = one;
            continue;
        }
        Matrix<T> trans(h + 1, h + 1, zero);
        for (std::size_t j = 0; j <= h; ++j) {
            for (std::size_t k = 0; k <= h; ++k) {
                if (R[i](j, k) == zero) continue;
                if (j == 0)
                    trans(j, k) = R[i](j, k);
                else
                    trans(j, k) = (one - exp(T(-lambda(i, j) * x[j - 1]))) * R[i](j, k);
                if (j != k) {
                    if (k == 0)
                        throw InputError(
                            "cache_ttl_lrua: the access graph routes into the not-cached node, "
                            "which has no characteristic time");
                    trans(k, j) = exp(T(-lambda(i, k) * x[k - 1]));
                }
            }
        }
        // Drop the nodes that no transition can reach (all-zero column).
        std::vector<std::size_t> chain;
        for (std::size_t c = 0; c <= h; ++c) {
            bool allzero = true;
            for (std::size_t r = 0; r <= h; ++r)
                if (!(trans(r, c) == zero)) {
                    allzero = false;
                    break;
                }
            if (!allzero) chain.push_back(c);
        }
        if (chain.empty()) throw NumericError("cache_ttl_lrua: the item has no reachable node");
        Matrix<T> P(chain.size(), chain.size());
        for (std::size_t a = 0; a < chain.size(); ++a)
            for (std::size_t b = 0; b < chain.size(); ++b) P(a, b) = trans(chain[a], chain[b]);
        const std::vector<T> ssprob = mc::dtmc_solve(P);

        std::vector<T> avgtime(chain.size());
        T den = zero;
        for (std::size_t a = 0; a < chain.size(); ++a) {
            const std::size_t node = chain[a];
            if (node == 0)
                avgtime[a] = one / lambda(i, 0);
            else
                avgtime[a] = (one - exp(T(-lambda(i, node) * x[node - 1]))) / lambda(i, node);
            den += ssprob[a] * avgtime[a];
        }
        if (den == zero) throw NumericError("cache_ttl_lrua: zero total holding time for an item");
        for (std::size_t a = 0; a < chain.size(); ++a)
            randprob(i, chain[a]) = ssprob[a] * avgtime[a] / den;
    }
    return randprob;
}

}  // namespace detail

/**
 * @param lambda (n x h+1) request rate of each item while it is at each node;
 *               column 0 is the rate while not cached
 * @param R      (n) access graphs, each (h+1 x h+1)
 * @param m      (h) list capacities
 * @param tol    relative tolerance on the characteristic times
 * @param maxswp cap on Gauss-Seidel sweeps
 * @return (n x h+1) time-stationary probabilities; column 0 is "not cached"
 */
template <class T>
Matrix<T> cache_ttl_lrua(const Matrix<T>& lambda, const std::vector<Matrix<T>>& R,
                         const std::vector<T>& m, const T& tol, unsigned maxswp = 200) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_ttl_lrua requires transcendental arithmetic");
    const std::size_t n = lambda.rows();
    const std::size_t h = m.size();
    if (n == 0) throw InputError("cache_ttl_lrua: no items");
    if (h == 0) throw InputError("cache_ttl_lrua: no lists");
    if (lambda.cols() != h + 1)
        throw InputError("cache_ttl_lrua: the rate matrix must have h+1 columns");
    if (R.size() != n) throw InputError("cache_ttl_lrua: one access graph per item is required");
    for (std::size_t i = 0; i < n; ++i)
        if (R[i].rows() != h + 1 || R[i].cols() != h + 1)
            throw InputError("cache_ttl_lrua: each access graph must be (h+1 x h+1)");
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t l = 0; l < h; ++l)
        if (!(m[l] > zero)) throw InputError("cache_ttl_lrua: list capacities must be positive");

    std::vector<T> x(h, num_traits<T>::from_int(1));
    for (unsigned sweep = 0; sweep < maxswp; ++sweep) {
        const std::vector<T> xold = x;
        for (std::size_t l = 0; l < h; ++l) {
            std::vector<T> work = x;
            auto resid = [&](const T& v) {
                work[l] = v;
                const Matrix<T> P = detail::lrua_randprob(lambda, R, work);
                T occ = zero;
                for (std::size_t i = 0; i < n; ++i) occ += P(i, l + 1);
                return T(occ - m[l]);
            };
            T lo = num_traits<T>::from_double(1e-12);
            T hi = x[l] > lo ? x[l] : num_traits<T>::from_int(1);
            if (!(resid(lo) < zero))
                throw NumericError("cache_ttl_lrua: list occupancy exceeds its capacity even at a "
                                   "vanishing characteristic time");
            // A LIST THAT CANNOT FILL HAS NO FINITE CHARACTERISTIC TIME, and
            // that is a state of the model rather than a numerical failure. The
            // occupancy is increasing in x and bounded by the number of items
            // the access graph can put in the list at a positive rate, so when
            // that bound is below the capacity the residual never changes sign
            // and the equation has no root: the list simply holds everything it
            // ever sees. It is REACHED IN PRACTICE by SolverLN, whose first
            // sweep solves the cache layer before any throughput has been
            // propagated into it -- every rate is then zero, the occupancy is
            // identically zero, and expanding the bracket only overflows.
            //
            // THE PROBE FOR THAT STATE MUST NOT UNDERFLOW, which is why it sits
            // at 700/max(lambda) rather than at any larger number. Once
            // exp(-lambda x) reaches exactly zero the item can no longer leave
            // list l, the embedded chain loses that edge, and dtmc_solve is
            // handed a chain whose recurrent class is a single absorbing node --
            // which ctmc_solve TRIMS, returning occupancy zero for the very list
            // that the limit fills completely. The probe would then read
            // "cannot fill" for a list that saturates, and the next sweep, run
            // with that value in hand, degenerates the chain outright. At
            // 700/max(lambda) every exponential is still a normal double
            // (exp(-700) ~ 1e-305), so the limit is evaluated to full precision
            // with the chain intact. The reference reaches the same point by a
            // different route: `cache_ttl_lrua.m` solves the same system with
            // FSOLVE, a least-squares method that needs no sign change and
            // returns its last iterate.
            T colmax = zero;
            for (std::size_t i = 0; i < n; ++i)
                if (lambda(i, l + 1) > colmax) colmax = lambda(i, l + 1);
            // No rate in the column: x[l] multiplies nothing, so nothing underflows.
            const T saturated = colmax > zero
                                    ? T(num_traits<T>::from_double(700.0) / colmax)
                                    : num_traits<T>::from_double(1e30);
            if (!(resid(saturated) > zero)) {
                x[l] = saturated;
                continue;
            }
            // Bracket by doubling, as bracket_expand does, but never past the
            // probe: the sign change is known to lie at or below it.
            const T two = num_traits<T>::from_int(2);
            if (!(hi < saturated)) hi = saturated;
            T fhi = resid(hi);
            while (fhi < zero) {
                hi *= two;
                if (!(hi < saturated)) {
                    hi = saturated;
                    break;
                }
                fhi = resid(hi);
            }
            const RootResult<T> r = root_bisect<T>(resid, lo, hi, T(tol * hi), 400);
            x[l] = r.root;
        }
        T rel = zero;
        for (std::size_t l = 0; l < h; ++l) {
            const T den = xold[l] > tol ? xold[l] : tol;
            const T d = num_abs(T(x[l] - xold[l])) / den;
            if (d > rel) rel = d;
        }
        if (rel < tol) break;
    }
    return detail::lrua_randprob(lambda, R, x);
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_TTL_LRUA_H
