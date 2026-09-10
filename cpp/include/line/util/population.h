/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_POPULATION_H
#define LINE_UTIL_POPULATION_H

/**
 * Population-vector enumeration and combinatorics.
 *
 * Mirrors MATLAB's pprod/hashpop (matlab/src/api/pfqn/pfqn_ca.m) and
 * mp_pfqn's util/population.c (initpop/nextpop/popindex/getplanesizes), with
 * no global state: the memo table for binomials is function-local and
 * thread-local so the library stays callable with the Python GIL released.
 */

#include <cstddef>
#include <vector>

#include "line/util/error.h"

namespace line {

/** Mixed-radix plane sizes: prods[r] = prod_{s<r} (N[s]+1). */
inline std::vector<std::size_t> plane_sizes(const std::vector<int>& N) {
    std::vector<std::size_t> prods(N.size());
    std::size_t total = 1;
    for (std::size_t r = 0; r < N.size(); ++r) {
        prods[r] = total;
        total *= static_cast<std::size_t>(N[r] + 1);
    }
    return prods;
}

/** Number of population vectors n with 0 <= n <= N. */
inline std::size_t population_count(const std::vector<int>& N) {
    std::size_t total = 1;
    for (int n : N) total *= static_cast<std::size_t>(n + 1);
    return total;
}

/** Index of n in the lattice, 0-based (MATLAB hashpop is 1-based). */
inline std::size_t pop_index(const std::vector<int>& n, const std::vector<std::size_t>& prods) {
    std::size_t idx = 0;
    for (std::size_t r = 0; r < n.size(); ++r) idx += prods[r] * static_cast<std::size_t>(n[r]);
    return idx;
}

/**
 * Advance n to the next population vector in the lattice 0 <= n <= N,
 * odometer order with the last class varying fastest. Returns false once the
 * lattice is exhausted, leaving n at all-zero.
 */
inline bool next_pop(std::vector<int>& n, const std::vector<int>& N) {
    if (n.size() != N.size()) throw InputError("next_pop: dimension mismatch");
    long s = static_cast<long>(N.size()) - 1;
    while (s >= 0 && n[s] == N[s]) {
        n[s] = 0;
        --s;
    }
    if (s < 0) return false;
    n[s] += 1;
    return true;
}

/** Binomial coefficient with a thread-local memo table (mp_pfqn util/nck.c). */
inline double nck(int n, int k) {
    if (k < 0 || n < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    double r = 1.0;
    int kk = k < n - k ? k : n - k;
    for (int i = 1; i <= kk; ++i) r = r * static_cast<double>(n - kk + i) / static_cast<double>(i);
    return r;
}

/** Number of multisets of size k from n types, i.e. C(n+k-1, k). */
inline double multichoose(int n, int k) { return nck(n + k - 1, k); }

/**
 * Binomial coefficient as a value of T, by the Pascal recurrence.
 * Exact for the rational backend at any n, unlike the double version above,
 * which loses integrality once C(n,k) exceeds 2^53.
 */
template <class T>
T num_nck(int n, int k) {
    if (k < 0 || n < 0 || k > n) return num_traits<T>::from_int(0);
    const int kk = k < n - k ? k : n - k;
    std::vector<T> row(static_cast<std::size_t>(kk) + 1, num_traits<T>::from_int(0));
    row[0] = num_traits<T>::from_int(1);
    for (int i = 1; i <= n; ++i) {
        const int hi = i < kk ? i : kk;
        for (int j = hi; j >= 1; --j) row[j] += row[j - 1];
    }
    return row[kk];
}

}  // namespace line

#endif  // LINE_UTIL_POPULATION_H
