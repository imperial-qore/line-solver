/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_PARAMS_H
#define LINE_API_MAPQN_MAPQN_PARAMS_H

/**
 * Model parameters and variable indexing shared by the mapqn QR bounds.
 *
 * A MAP queueing network here is M queues, N circulating jobs, K(i) phases at
 * queue i, completion rates mu{i}(k,h), background (non-completion) phase
 * transition rates v{i}(k,h), load-dependent scalings alpha(i,n) and routing
 * probabilities r(i,j). This mirrors the params struct that
 * matlab/lib/qrf/mapqn_bnd_qr_ld.m takes and the Mapqn_parameters class in
 * jar/src/main/java/jline/api/mapqn/.
 *
 * Index convention: queues and phases are 0-based here (C++ convention),
 * populations are as written, 0..N. The MATLAB reference is 1-based in queues
 * and phases, so every loop `for i = 1:M` becomes `for i = 0; i < M; ++i` and
 * every `q_val(i,j,k,h,n)` argument drops by one except n. This matches the
 * JAR, whose q takes 0-based i,j,k,h and a 1-based population n (see the mapqn
 * section of _kb/03-api-layer.md).
 *
 * Arithmetic: nothing here is transcendental. q is a sum of products of model
 * data, so at T = line::Rational every LP coefficient is an exact rational and
 * the bound the LP returns is the exact optimum of the exact polytope.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mapqn {

/** Parameters of a MAP queueing network for the QR bounds. */
template <class T>
struct MapqnParams {
    int M = 0;                    ///< number of queues
    int N = 0;                    ///< total population
    std::vector<int> K;           ///< K[i] = number of phases at queue i
    std::vector<Matrix<T>> mu;    ///< mu[i] is K(i) x K(i), completion rates
    std::vector<Matrix<T>> v;     ///< v[i] is K(i) x K(i), background rates
    Matrix<T> alpha;              ///< M x N load-dependent scalings, 1 beyond the table
    Matrix<T> r;                  ///< M x M routing probabilities
    T Z = T();                    ///< think time (delay model only)
    T D1 = num_traits<T>::from_int(1);  ///< service demand at queue 1 (delay model only)

    void validate() const {
        if (M <= 0) throw InputError("mapqn: M must be positive");
        if (N < 0) throw InputError("mapqn: N must be nonnegative");
        if (static_cast<int>(K.size()) != M) throw InputError("mapqn: K has the wrong length");
        if (static_cast<int>(mu.size()) != M || static_cast<int>(v.size()) != M)
            throw InputError("mapqn: mu and v must have one entry per queue");
        for (int i = 0; i < M; ++i) {
            if (K[i] <= 0) throw InputError("mapqn: every queue needs at least one phase");
            const std::size_t k = static_cast<std::size_t>(K[i]);
            if (mu[i].rows() != k || mu[i].cols() != k)
                throw InputError("mapqn: mu{i} must be K(i) x K(i)");
            if (v[i].rows() != k || v[i].cols() != k)
                throw InputError("mapqn: v{i} must be K(i) x K(i)");
        }
        if (r.rows() != static_cast<std::size_t>(M) || r.cols() != static_cast<std::size_t>(M))
            throw InputError("mapqn: r must be M x M");
        if (!alpha.empty() && alpha.rows() != static_cast<std::size_t>(M))
            throw InputError("mapqn: alpha must have M rows");
    }
};

/**
 * q(i,j,k,h,n): rate at which queue i, holding n jobs and in phase k, moves to
 * phase h while routing a job to queue j.
 *
 * Port of the local q_func in mapqn_bnd_qr_ld.m. n == 0 returns 0 (an empty
 * queue completes nothing); the load-dependent factor alpha(i,n) is 1 beyond
 * the width of the alpha table, exactly as MATLAB's `if n <= size(alpha,2)`
 * guard does. The i == j branch adds the background rate v, because a
 * self-routing completion and a phase change without completion are
 * indistinguishable in the marginal process.
 */
template <class T>
T mapqn_q(const MapqnParams<T>& p, int i, int j, int k, int h, int n) {
    if (n == 0) return T();
    T alpha_val = num_traits<T>::from_int(1);
    if (!p.alpha.empty() && static_cast<std::size_t>(n) <= p.alpha.cols())
        alpha_val = p.alpha(static_cast<std::size_t>(i), static_cast<std::size_t>(n - 1));
    const std::size_t ki = static_cast<std::size_t>(k), hi = static_cast<std::size_t>(h);
    if (j != i) {
        const T base = p.r(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) * p.mu[i](ki, hi);
        return T(base * alpha_val);
    }
    const T base = p.v[i](ki, hi) +
                   p.r(static_cast<std::size_t>(i), static_cast<std::size_t>(i)) * p.mu[i](ki, hi);
    return T(base * alpha_val);
}

/**
 * Flat index of the joint variable p2(j,nj,k,i,ni,h).
 *
 * The enumeration order is the one mapqn_bnd_qr_ld.m builds:
 *   for j, for nj = 0..N, for k = 1..K(j), for i, for ni = 0..N, for h = 1..K(i)
 * which factorizes as outer(j,nj,k) * B + inner(i,ni,h) with the same map used
 * for both halves and B = (N+1) * sum_i K(i). Keeping the two halves identical
 * is what makes the SYMMETRY family a plain index swap.
 */
struct P2Index {
    int M = 0, N = 0;
    std::vector<int> K;
    std::vector<int> cumK;  ///< cumK[i] = sum_{i' < i} K(i')
    std::size_t block = 0;  ///< (N+1) * sum_i K(i)

    P2Index() {}
    P2Index(int m, int n, const std::vector<int>& k) : M(m), N(n), K(k) {
        cumK.assign(static_cast<std::size_t>(M) + 1, 0);
        for (int i = 0; i < M; ++i) cumK[i + 1] = cumK[i] + K[i];
        block = static_cast<std::size_t>(N + 1) * static_cast<std::size_t>(cumK[M]);
    }

    /** Half-index of (queue, population, phase); identical for both halves. */
    std::size_t half(int i, int ni, int h) const {
        return static_cast<std::size_t>(N + 1) * static_cast<std::size_t>(cumK[i]) +
               static_cast<std::size_t>(ni) * static_cast<std::size_t>(K[i]) +
               static_cast<std::size_t>(h);
    }

    std::size_t operator()(int j, int nj, int k, int i, int ni, int h) const {
        return half(j, nj, k) * block + half(i, ni, h);
    }

    std::size_t num_vars() const { return block * block; }
};

/** Result of a QR bound solve. */
template <class T>
struct MapqnQrResult {
    bool ok = false;             ///< the LP reached an optimal vertex
    std::string status;          ///< textual LP status
    T objective = T();           ///< the bound
    std::vector<T> x;            ///< full solution vector, indexed by P2Index
    std::vector<Matrix<T>> p2marginals;  ///< p2marginals[j](nj, k) = p2(j,nj,k,j,nj,k)
    std::size_t num_vars = 0;
    std::size_t num_rows = 0;
    std::size_t iterations = 0;
};

/** Which direction the bound is taken in. */
enum class MapqnSense { Max, Min };

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_PARAMS_H
