/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_TYPES_H
#define LINE_API_MDD_MDD_TYPES_H

/**
 * The rate side of the decision-diagram domain: local matrices, events, the
 * Kronecker descriptor, and the options/result of the level aggregation.
 *
 * Port of the JAR classes MddLocalMatrix, MddEvent, MddServiceLaw,
 * MddDescriptor, MddMcdOptions and MddMcdResult, which the MATLAB and python
 * twins carry as struct fields of the same names.
 *
 * The rate matrix of a structured model is R = sum_e (kron_k W_k^e) restricted
 * to the reachable set, with W_k^e[i,j] = lambda_k^e[i] * Prob_k^e(i,j) (Eq. 1
 * of Miner-Ciardo-Donatelli, SIGMETRICS 2000). An event touches only the levels
 * it names; every other level carries the identity, which `mdd_mcd` supplies
 * rather than storing.
 */

#include <cstddef>
#include <functional>
#include <map>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mdd {

/**
 * A local rate matrix W_k^e of the Kronecker descriptor, held row-compressed.
 *
 * The aggregation only ever walks a row of W and needs its row sums (the local
 * enabling rates lambda_k^e), so the rows are stored as parallel column/value
 * arrays rather than as a general sparse matrix. Duplicate triplets are summed
 * when the matrix is built, so a caller may emit the same (i,j) more than once.
 */
template <class T>
struct MddLocalMatrix {
    /** Order of the (square) local matrix, i.e. the level domain. */
    std::size_t dim = 0;
    /** cols[i] holds the column indices of the nonzeros of row i. */
    std::vector<std::vector<std::size_t>> cols;
    /** vals[i] holds the values of the nonzeros of row i, aligned with cols[i]. */
    std::vector<std::vector<T>> vals;
    /** row_sum[i] is the local enabling rate lambda[i]. */
    std::vector<T> row_sum;
    /** Total number of stored nonzeros. */
    std::size_t nnz = 0;

    /** The identity of the given order, used for a level an event does not touch. */
    static MddLocalMatrix<T> identity(std::size_t d) {
        const T one = num_traits<T>::from_int(1);
        MddLocalMatrix<T> m;
        m.dim = d;
        m.cols.assign(d, std::vector<std::size_t>());
        m.vals.assign(d, std::vector<T>());
        m.row_sum.assign(d, one);
        for (std::size_t i = 0; i < d; ++i) {
            m.cols[i].push_back(i);
            m.vals[i].push_back(one);
        }
        m.nnz = d;
        return m;
    }

    /** Incremental triplet builder; duplicate entries are accumulated. */
    class Builder {
    public:
        explicit Builder(std::size_t d) : dim_(d), rows_(d) {}

        /** Accumulate value into entry (i,j). A zero value is dropped. */
        Builder& add(std::size_t i, std::size_t j, const T& value) {
            const T zero = num_traits<T>::from_int(0);
            if (value == zero) return *this;
            if (i >= dim_ || j >= dim_)
                throw InputError("MddLocalMatrix::Builder: index outside the level domain");
            typename std::map<std::size_t, T>::iterator it = rows_[i].find(j);
            if (it == rows_[i].end())
                rows_[i][j] = value;
            else
                it->second = T(it->second + value);
            return *this;
        }

        MddLocalMatrix<T> build() const {
            const T zero = num_traits<T>::from_int(0);
            MddLocalMatrix<T> m;
            m.dim = dim_;
            m.cols.assign(dim_, std::vector<std::size_t>());
            m.vals.assign(dim_, std::vector<T>());
            m.row_sum.assign(dim_, zero);
            for (std::size_t i = 0; i < dim_; ++i) {
                T s = zero;
                // std::map iterates in ascending key order, so two builds of the
                // same matrix agree entry for entry
                for (typename std::map<std::size_t, T>::const_iterator it = rows_[i].begin();
                     it != rows_[i].end(); ++it) {
                    m.cols[i].push_back(it->first);
                    m.vals[i].push_back(it->second);
                    s += it->second;
                }
                m.row_sum[i] = s;
                m.nnz += rows_[i].size();
            }
            return m;
        }

    private:
        std::size_t dim_;
        std::vector<std::map<std::size_t, T>> rows_;
    };
};

/** One event of the Kronecker rate descriptor. */
template <class T>
struct MddEvent {
    /** Station (or transition node) the event departs from, 0-based. */
    std::size_t a = 0;
    /** Station (or mode) the event arrives at, 0-based; equals a for an internal event. */
    std::size_t b = 0;
    /** Levels the event touches, as 0-based level indices, aligned with W. */
    std::vector<std::size_t> lev;
    /** Local matrices at the levels named by lev. */
    std::vector<MddLocalMatrix<T>> W;
};

/**
 * Phase-type service law of one station, as a Markovian (D0,D1) pair.
 *
 * D0 holds the phase transitions that do not complete a service and D1 those
 * that do, so the exit-rate vector is t0 = D1*1. For a renewal law D1 = t0*pie,
 * and the entry law pie is then DERIVED from D1 rather than assumed (see
 * `mdd_entry_law`). Set `pie` only to override that.
 */
template <class T>
struct MddServiceLaw {
    Matrix<T> D0;
    Matrix<T> D1;
    /** Optional entry law; empty to derive it from D1. */
    std::vector<T> pie;
    bool present = false;  ///< false marks "this station is exponential"

    std::size_t phases() const { return D0.rows(); }
};

/** Successor function over local indices, for `mdd_reachset`. */
typedef std::function<std::vector<std::vector<int>>(const std::vector<int>&)> MddNextState;

/**
 * Kronecker rate descriptor of a structured model, the input of `mdd_mcd`.
 *
 * Built by `mdd_descriptor` (count-plus-in-service-phase local states,
 * non-preemptive), `mdd_ps` (per-phase-count local states, shared servers) or
 * `spn::spn_mdd` (a stochastic Petri net).
 */
template <class T>
struct MddDescriptor {
    /** Number of levels, i.e. stations or places. */
    std::size_t K = 0;
    /** Closed population; the conservation law the level marginals must satisfy. */
    int N = 0;
    /** Local domain per level. */
    std::vector<int> domain;
    /** Station service rates, 1/E[S]; empty for a descriptor with no queueing parameters. */
    std::vector<T> mu;
    /** Servers per station; infinite for a delay station. */
    std::vector<double> servers;
    /** Station-to-station routing matrix. */
    Matrix<T> P;
    /** Phases per station, 1 when exponential. */
    std::vector<std::size_t> nphases;
    /**
     * valuemap[i][idx] is the physical occupancy of level i in local state idx.
     *
     * A level whose local state encodes more than a count (a station holding
     * both a population and a service phase) needs this map; without one the
     * index would be the quantity.
     */
    std::vector<std::vector<double>> valuemap;
    /** Initial local index per level. */
    std::vector<int> init;
    /** Successor function over local indices. */
    MddNextState nextfun;
    /** The events of the descriptor. */
    std::vector<MddEvent<T>> events;
    /**
     * Optional conservation law as weights' * QLen = value, overriding the
     * closed-population test. Empty when the population N is the invariant.
     */
    std::vector<double> invariant_weights;
    /** Value of the invariant when invariant_weights is set. */
    double invariant_value = 0.0;
};

/**
 * Knobs of the level iteration in `mdd_mcd`.
 *
 * The defaults are deliberately much tighter than a solver-level fixed-point
 * tolerance: the level iteration is an INNER numerical solve and `mdd_mcd`
 * verifies the population invariant at 1e-6, so a loose tolerance converges
 * short of the fixed point and trips that guard. Do not wire an AMVA-sized
 * iter_tol into these.
 */
struct MddMcdOptions {
    /** Convergence tolerance on the level marginals. */
    double tol = 1e-12;
    /** Maximum coupled sweeps before the iteration is declared non-convergent. */
    int maxiter = 500;
    /**
     * The reference's 'verbose' knob is NOT carried: it is a console trace of
     * the level sizes and the iteration count, and this port keeps the api layer
     * silent, as `infer_lqn_ekf` does. `MddMcdResult` returns level_sizes and
     * iters, so the same numbers are available to a caller that wants them.
     */
    /** Optional warm-start level vectors, one per paper level; empty to start uniform. */
    std::vector<std::vector<double>> initpik;
};

/** Result of the Miner-Ciardo-Donatelli level aggregation. */
template <class T>
struct MddMcdResult {
    /** Mean occupancy per station (or place), in station order. */
    std::vector<T> QLen;
    /** Per-station throughput; empty when the descriptor carries no queueing parameters. */
    std::vector<T> X;
    /** Per-station utilization; empty when the descriptor carries no queueing parameters. */
    std::vector<T> U;
    /** pik[k] is the level-k stationary vector over M_k, in paper orientation. */
    std::vector<std::vector<T>> pik;
    /** Mrows[k][r] = {node id, local value} of row r of M_k. */
    std::vector<std::vector<std::pair<int, int>>> Mrows;
    /** |M_k| per paper level. */
    std::vector<std::size_t> level_sizes;
    /** Fixed-point iterations performed. */
    int iters = 0;
    /**
     * max |A(p)| per paper level: the largest number of distinct root-to-node
     * paths at that level. 1 means no node there is shared, so conditioning on
     * the node equals conditioning on the whole path above it.
     */
    std::vector<double> paths_per_level;
    /**
     * True certifies the result is EXACT with no reference solve needed; false
     * means "not certified by this test", never "approximate" -- a product-form
     * model is exact however much its diagram shares.
     */
    bool no_aggregation = false;
};

/**
 * Entry law of a phase-type station, taken as given or derived from D1.
 *
 * A {D0,D1} pair carries its own restart law: for a renewal process D1 = t0*pie,
 * so every row with a positive exit rate is proportional to pie. Deriving it is
 * not optional -- defaulting to e_1 instead silently replaces a hyperexponential
 * (whose D0 is diagonal, so a job entering phase 1 can never leave it) by an
 * exponential at the phase-1 rate.
 */
template <class T>
std::vector<T> mdd_entry_law(const std::vector<T>& given, const Matrix<T>& D1, std::size_t h,
                             std::size_t i, const std::string& caller) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (!given.empty()) {
        std::vector<T> v = given;
        T s = zero;
        for (std::size_t a = 0; a < v.size(); ++a) s += v[a];
        for (std::size_t a = 0; a < v.size(); ++a) v[a] = T(v[a] / s);
        return v;
    }
    std::vector<T> t0(h, zero);
    long live = -1;
    std::size_t nlive = 0;
    for (std::size_t a = 0; a < h; ++a) {
        T s = zero;
        for (std::size_t b = 0; b < h; ++b) s += D1(a, b);
        t0[a] = s;
        if (s > zero) {
            ++nlive;
            if (live < 0) live = static_cast<long>(a);
        }
    }
    if (live < 0) {
        std::vector<T> pie(h, zero);
        pie[0] = one;
        return pie;
    }
    const std::size_t l = static_cast<std::size_t>(live);
    std::vector<T> first(h, zero);
    for (std::size_t b = 0; b < h; ++b) first[b] = T(D1(l, b) / t0[l]);
    // A non-renewal MAP restarts in a law that depends on the phase it left
    // from, which one entry vector cannot express; the composite level would
    // silently model the renewal process instead.
    if (nlive > 1) {
        for (std::size_t a = l + 1; a < h; ++a) {
            if (!(t0[a] > zero)) continue;
            for (std::size_t b = 0; b < h; ++b) {
                const double d = num_traits<T>::to_double(T(D1(a, b) / t0[a])) -
                                 num_traits<T>::to_double(first[b]);
                if (d > 1e-9 || d < -1e-9)
                    throw InputError(caller + ": station " + std::to_string(i + 1) +
                                     " carries a service law whose restart distribution depends "
                                     "on the completing phase (a non-renewal MAP); the local "
                                     "state names one entry law, so this encoding cannot "
                                     "represent it");
            }
        }
    }
    return first;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_TYPES_H
