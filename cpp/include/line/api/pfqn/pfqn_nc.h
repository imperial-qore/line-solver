/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_NC_H
#define LINE_API_PFQN_NC_H

/**
 * Normalizing constant of a product-form queueing network: the dispatcher.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_nc.m. The routine is two things
 * at once, and they are kept separate here:
 *
 *  1. a MODEL REDUCTION, which folds the open classes into rescaled demands,
 *     drops empty and demand-free classes and stations, normalizes the demands
 *     per class, and closes the degenerate cases (single station, M identical
 *     replicas, delay only) in closed form. This is where most of the file is,
 *     and it is exact.
 *  2. a DISPATCH to one of the normalizing-constant algorithms.
 *
 * OPEN CLASSES. MATLAB marks an open class by N_r = Inf. There is no infinity
 * in an exact field, so this port marks it by a NEGATIVE population, which is
 * also what the .qn interchange format uses. Each station's demands are
 * inflated by 1/(1 - sum_r lambda_r L(i,r)) and the open queue lengths are
 * read off directly, exactly as in the reference. Note that the reference
 * leaves lGopen at zero (the line accumulating sum log Ut is commented out),
 * so the returned constant is the CLOSED-CONDITIONAL one; that convention is
 * preserved, because changing it would silently rescale every caller's result.
 *
 * ARITHMETIC, NOT ALGORITHMS, IS WHAT IS REFUSED. The whole ladder of
 * `compute_norm_const` is dispatched here: ca, clw, cub/gm, kt, bkt, lekt, le, ble, dir, aghq, ls, is,
 * mci, imci, sampling, mmint2/gleint, pana, propfair, rgf, mva, exact, comom
 * and recal. Every one of them but ca / exact / recal / mva / comom is an
 * asymptotic or Monte Carlo estimator formed in logarithms, so it is compiled
 * only when `num_traits<T>::has_transcendental` and refused BY NAME otherwise
 * rather than silently redirected to the convolution: a normalizing constant
 * computed by an algorithm the caller did not ask for is indistinguishable
 * from the right answer until it is wrong. NcMethod::Default's multi-station
 * branch (cub for sum(N) < 1e3, le above it) is refused on the same grounds in
 * exact arithmetic. `rgf` is the one two-sided name: on a SINGLE-CLASS model it
 * is the log-domain generating-function recursion and is refused in an exact
 * field, while a multiclass request reduces to ca and stays exact, reported as
 * 'rgf/ca'.
 *
 * Arithmetic: EXACT-CAPABLE for methods Ca, Exact, Recal, Mva, Comom and for
 * the whole reduction. The reference's tolerance-based filters
 * (options.tol, GlobalConstants.FineTol) become exact zero tests, for the same
 * reason as in pfqn_unique: a tolerant filter perturbs the constant and has no
 * meaning in the rational field. Pass a positive atol to recover the tolerant
 * behaviour in any arithmetic.
 *
 * WHY THE INEXACT BRANCHES REPORT lG AND NOT G. An estimator returns the LOG
 * constant; exponentiating it and taking the log again loses the answer
 * outright once lG passes ~709. Those branches therefore accumulate lG
 * additively (log Gscale + log Gzdem + lG_core) and derive G from it, while
 * the exact branches keep the exact product and take its log, as before.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_clw.h"
#include "line/api/pfqn/pfqn_comomrm.h"
#include "line/api/pfqn/pfqn_cub.h"
#include "line/api/pfqn/pfqn_cub_evals.h"
#include "line/api/pfqn/pfqn_gerasimov.h"
#include "line/api/pfqn/pfqn_is.h"
#include "line/api/pfqn/pfqn_bk.h"
#include "line/api/pfqn/pfqn_kt.h"
#include "line/api/pfqn/pfqn_bkt.h"
#include "line/api/pfqn/pfqn_lekt.h"
#include "line/api/pfqn/pfqn_le.h"
#include "line/api/pfqn/pfqn_aghq.h"
#include "line/api/pfqn/pfqn_ble.h"
#include "line/api/pfqn/pfqn_ls.h"
#include "line/api/pfqn/pfqn_mci.h"
#include "line/api/pfqn/pfqn_mcmc.h"
#include "line/api/pfqn/pfqn_mmint2.h"
#include "line/api/pfqn/pfqn_mmsample2.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_panacea.h"
#include "line/api/pfqn/pfqn_propfair.h"
#include "line/api/pfqn/pfqn_recal.h"
#include "line/api/pfqn/pfqn_explicit.h"
#include "line/api/pfqn/pfqn_rgf.h"
#include "line/api/pfqn/pfqn_rgfmc.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** The methods this port dispatches, one per `compute_norm_const` case. */
enum class NcMethod {
    Default,
    Adaptive,  ///< the reference groups 'adaptive' with 'default'
    Ca,
    Exact,
    Recal,
    Mva,
    Comom,
    Clw,
    Cub,
    Gm,  ///< the reference's alias of 'cub'
    Kt,
    Bkt,      ///< KT minus the exact Stirling remainder of each Laplaced class (BKT)
    Lekt,      ///< the estimator Ble and Bkt both compute, on the cheaper side
    Bk,        ///< Birman-Kogan saddle point with bottleneck detection
    Bkue,      ///< Birman-Kogan uniform (van der Waerden) expansion, single chain
    Lc,        ///< Birman-Kogan Algorithm 2, single chain subproblems by MVA
    LcUe,      ///< Algorithm 2 with the uniform expansion as the single chain solver
    Le,
    Ble,  ///< LE plus the empirical eps->0 correction
    Aghq,  ///< adaptive Gauss-Hermite over the simplex; q=1 is Le
    Ls,
    Is,
    Mci,
    Imci,
    Mcmc,  ///< Chen-O'Cinneide regularization; supplies X and Q, never a constant
    Sampling,
    Mmint2,
    Gleint,  ///< the reference's alias of 'mmint2'
    Pana,
    Propfair,
    Rgf,  ///< recursion by generating functions; residues beyond one class
    Divdiff,  ///< divided-difference closed form; no think time, no load dependence
    Ger  ///< residue closed form; free in the eliminated class populations
};

inline const char* nc_method_name(NcMethod m) {
    switch (m) {
        case NcMethod::Default: return "default";
        case NcMethod::Adaptive: return "adaptive";
        case NcMethod::Ca: return "ca";
        case NcMethod::Exact: return "exact";
        case NcMethod::Recal: return "recal";
        case NcMethod::Mva: return "mva";
        case NcMethod::Comom: return "comom";
        case NcMethod::Clw: return "clw";
        case NcMethod::Cub: return "cub";
        case NcMethod::Gm: return "gm";
        case NcMethod::Kt: return "kt";
        case NcMethod::Bkt: return "bkt";
        case NcMethod::Lekt: return "lekt";
        case NcMethod::Bk: return "bk";
        case NcMethod::Bkue: return "bkue";
        case NcMethod::Lc: return "lc";
        case NcMethod::LcUe: return "lc.ue";
        case NcMethod::Le: return "le";
        case NcMethod::Ble: return "ble";
        case NcMethod::Aghq: return "aghq";
        case NcMethod::Ls: return "ls";
        case NcMethod::Is: return "is";
        case NcMethod::Mci: return "mci";
        case NcMethod::Imci: return "imci";
        case NcMethod::Mcmc: return "mcmc";
        case NcMethod::Sampling: return "sampling";
        case NcMethod::Mmint2: return "mmint2";
        case NcMethod::Gleint: return "gleint";
        case NcMethod::Pana: return "pana";
        case NcMethod::Propfair: return "propfair";
        case NcMethod::Rgf: return "rgf";
        case NcMethod::Divdiff: return "divdiff";
        case NcMethod::Ger: return "ger";
    }
    return "default";
}

/** Map a method name to its enum; throws UnsupportedError on an unknown one. */
inline NcMethod nc_method_of(const std::string& s) {
    if (s == "default") return NcMethod::Default;
    if (s == "adaptive") return NcMethod::Adaptive;
    if (s == "ca") return NcMethod::Ca;
    if (s == "exact") return NcMethod::Exact;
    if (s == "recal") return NcMethod::Recal;
    if (s == "mva") return NcMethod::Mva;
    if (s == "comom") return NcMethod::Comom;
    if (s == "clw") return NcMethod::Clw;
    if (s == "cub") return NcMethod::Cub;
    if (s == "gm") return NcMethod::Gm;
    if (s == "kt") return NcMethod::Kt;
    if (s == "bkt") return NcMethod::Bkt;
    if (s == "lekt") return NcMethod::Lekt;
    if (s == "bk") return NcMethod::Bk;
    if (s == "bkue") return NcMethod::Bkue;
    if (s == "lc") return NcMethod::Lc;
    if (s == "lc.ue") return NcMethod::LcUe;
    if (s == "le") return NcMethod::Le;
    if (s == "ble") return NcMethod::Ble;
    if (s == "aghq") return NcMethod::Aghq;
    if (s == "ls") return NcMethod::Ls;
    if (s == "is") return NcMethod::Is;
    if (s == "mci") return NcMethod::Mci;
    if (s == "imci") return NcMethod::Imci;
    if (s == "mcmc") return NcMethod::Mcmc;
    if (s == "sampling") return NcMethod::Sampling;
    if (s == "mmint2") return NcMethod::Mmint2;
    if (s == "gleint") return NcMethod::Gleint;
    if (s == "pana") return NcMethod::Pana;
    if (s == "propfair") return NcMethod::Propfair;
    if (s == "rgf") return NcMethod::Rgf;
    if (s == "divdiff") return NcMethod::Divdiff;
    if (s == "ger") return NcMethod::Ger;
    throw UnsupportedError("pfqn_nc: unrecognized method '" + s + "'");
}

/** The `options` fields `compute_norm_const` reads beyond the method itself. */
struct NcOptions {
    std::size_t samples = 100000;  ///< SolverOptions('NC').samples
    unsigned long seed = 23000;    ///< SolverOptions('NC').seed
    double tol = 1e-6;             ///< handed to pfqn_comomrm
    /// options.config.aghq_nodes: nodes per simplex direction of the adaptive
    /// Gauss-Hermite rule. q = 1 reproduces pfqn_le and the rule costs q^(M-1)
    /// evaluations, so the default stays small.
    std::size_t aghq_nodes = 3;
    /// options.config.mcmc_batches: batches pfqn_mcmc splits its run into
    std::size_t mcmc_batches = MCMC_DEFAULT_BATCHES;
    /// options.config.mcmc_burnin: warm-up fraction pfqn_mcmc discards
    double mcmc_burnin = MCMC_DEFAULT_BURNIN;
};

/**
 * Refuse a method in an arithmetic it has no meaning in.
 * Kept as a function so the message is identical wherever it is raised.
 */
inline void pfqn_nc_refuse(const std::string& method) {
    throw UnsupportedError("pfqn_nc: method '" + method +
                           "' is an asymptotic or Monte Carlo estimator formed in logarithms and "
                           "needs transcendental arithmetic. "
                           "Use 'ca', 'exact', 'recal', 'mva' or 'comom' for an exact constant.");
}

template <class T>
struct NcDispatchResult {
    T G;                ///< normalizing constant, exact when the method is
    double lG;          ///< its logarithm
    std::vector<T> X;   ///< (R) per-class throughput, empty when not produced
    Matrix<T> Q;        ///< (M x R) queue lengths, empty when not produced
    std::string method; ///< the algorithm actually used
    /**
     * False when the METHOD DECLINED THE MODEL, which is the reference's
     * `lG = []`: pana outside normal usage, or mmint2/gleint on a model with
     * more than one queueing station. MATLAB warns and returns empty, and
     * `getAvg` then renders a table of ZEROS while still reporting a completed
     * analysis; the caller reproduces that by returning an all-zero solution.
     *
     * This is NOT the same as a refusal. A refusal throws, and the reference
     * throws too -- `comom` on a multi-queue model is `line_error` in
     * `pfqn_nc.m`, added deliberately with the comment that "a silently zeroed
     * result is worse than no result". Only the two branches above decline.
     * Ruled by the user on 2026-07-25 (register row N1).
     */
    bool valid = true;
};

/**
 * @param lambda (R) arrival rates; zero on the closed classes, may be empty
 * @param L      (M x R) service demands
 * @param N      (R) populations; a NEGATIVE entry marks an open class
 * @param Z      (K x R) think times, summed over rows
 * @param method requested algorithm
 * @param atol   threshold below which a demand counts as zero; 0 for exact
 * @param nopt   sample count, seed and tolerance the estimators read
 */
template <class T>
NcDispatchResult<T> pfqn_nc(const std::vector<T>& lambda, const Matrix<T>& L,
                            const std::vector<int>& N, const Matrix<T>& Z, NcMethod method,
                            const T& atol, const NcOptions& nopt) {
    const std::size_t R = N.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    NcDispatchResult<T> res;
    res.G = one;
    res.lG = 0.0;
    // Every return BEFORE the dispatch reports 'exact', because every one of
    // them is a closed form rather than an algorithm: the reference sets this
    // up front and lets the dispatch overwrite it.
    res.method = "exact";

    if (R == 0) {
        res.G = zero;
        res.lG = -std::numeric_limits<double>::infinity();
        return res;
    }
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_nc: L and N disagree on the class count");

    const std::size_t M0 = L.empty() ? 0 : L.rows();
    std::vector<T> lam(R, zero);
    for (std::size_t r = 0; r < R && r < lambda.size(); ++r) lam[r] = lambda[r];

    // ---- open classes: inflate the demands, read off the open queue lengths --
    Matrix<T> Lw = L;
    Matrix<T> Qopen(M0, R, zero);
    std::vector<std::size_t> ocl;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] < 0) ocl.push_back(r);
    for (std::size_t i = 0; i < M0; ++i) {
        T u = one;
        for (std::size_t r = 0; r < R; ++r) u -= lam[r] * L(i, r);
        if (u == zero) {
            // A station the open classes fill exactly. This is not a pathology:
            // it is what a SOURCE row looks like, whose chain demand is 1/lambda
            // by construction. The reference divides by the zero, gets Inf and
            // NaN, and its demand filter (Lmax./Lsum > FineTol, false on NaN)
            // drops the row a few lines later. Dropping it here reaches the same
            // constant without carrying non-finite values through the reduction.
            for (std::size_t r = 0; r < R; ++r) Lw(i, r) = zero;
            continue;
        }
        for (std::size_t r = 0; r < R; ++r) {
            Lw(i, r) = L(i, r) / u;
            Qopen(i, r) = lam[r] * Lw(i, r) / u;
        }
    }

    // Closed populations, with the open classes zeroed out.
    std::vector<int> Nc(R, 0);
    for (std::size_t r = 0; r < R; ++r) Nc[r] = N[r] > 0 ? N[r] : 0;
    long Ntot = 0;
    for (int v : Nc) Ntot += v;
    if (Ntot == 0) return res;  // lG = 0, G = 1

    const auto colsum = [&](const Matrix<T>& A, std::size_t r) {
        T s = zero;
        for (std::size_t i = 0; i < A.rows(); ++i) s += A(i, r);
        return s;
    };

    // ---- drop the empty classes ---------------------------------------------
    std::vector<std::size_t> nnz;
    for (std::size_t r = 0; r < R; ++r)
        if (Nc[r] > 0) nnz.push_back(r);
    const std::size_t R1 = nnz.size();

    // per-class rescaling-order rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<T> scalevec(R1, one);
    Matrix<T> L1(M0, R1), Z1(Z.empty() ? 0 : Z.rows(), R1);
    for (std::size_t k = 0; k < R1; ++k) {
        const std::size_t r = nnz[k];
        T mx = zero;
        for (std::size_t i = 0; i < M0; ++i)
            if (Lw(i, r) > mx) mx = Lw(i, r);
        for (std::size_t i = 0; i < Z1.rows(); ++i)
            if (Z(i, r) > mx) mx = Z(i, r);
        if (mx > zero) scalevec[k] = mx;
        for (std::size_t i = 0; i < M0; ++i) L1(i, k) = Lw(i, r) / scalevec[k];
        for (std::size_t i = 0; i < Z1.rows(); ++i) Z1(i, k) = Z(i, r) / scalevec[k];
    }
    // Gscale is the factor the rescaling removed: prod_r scalevec_r^{N_r}.
    T Gscale = one;
    for (std::size_t k = 0; k < R1; ++k)
        Gscale *= num_pow_int(scalevec[k], static_cast<unsigned>(Nc[nnz[k]]));
    // The LOG path must never take the log of that product: scalevec = 3 with
    // N = 600 leaves double range and turned a perfectly good lG = 862 into inf.
    // Accumulate sum_r N_r log(scalevec_r) directly instead. The exact branches
    // keep using Gscale itself, where the product is representable by
    // construction.
    double lGscale = 0.0;
    if (num_traits<T>::has_transcendental) {
        for (std::size_t k = 0; k < R1; ++k)
            lGscale += static_cast<double>(Nc[nnz[k]]) *
                       std::log(num_traits<T>::to_double(scalevec[k]));
    }

    // ---- drop the stations with no demand at all -----------------------------
    std::vector<std::size_t> demSt;
    for (std::size_t i = 0; i < M0; ++i) {
        T rs = zero;
        for (std::size_t k = 0; k < R1; ++k) rs += L1(i, k);
        if (rs > atol) demSt.push_back(i);
    }
    Matrix<T> L2(demSt.size(), R1);
    for (std::size_t a = 0; a < demSt.size(); ++a)
        for (std::size_t k = 0; k < R1; ++k) L2(a, k) = L1(demSt[a], k);
    const std::size_t M = demSt.size();

    std::vector<int> N2(R1, 0);
    for (std::size_t k = 0; k < R1; ++k) N2[k] = Nc[nnz[k]];

    // Aggregate think times, one row.
    Matrix<T> Zag(Z1.rows() == 0 ? 0 : 1, R1);
    T Ztot = zero;
    for (std::size_t k = 0; k < R1; ++k) {
        T s = zero;
        for (std::size_t i = 0; i < Z1.rows(); ++i) s += Z1(i, k);
        if (Zag.rows() > 0) Zag(0, k) = s;
        Ztot += s;
    }

    // Delay-only constant prod_r Z_r^{N_r}/N_r!, needed by several branches.
    const auto delayG = [&](const std::vector<std::size_t>& cls) {
        T g = one;
        for (std::size_t k : cls) {
            T zs = zero;
            for (std::size_t i = 0; i < Z1.rows(); ++i) zs += Z1(i, k);
            g *= num_pow_int(zs, static_cast<unsigned>(N2[k])) /
                 num_factorial<T>(static_cast<unsigned>(N2[k]));
        }
        return g;
    };
    // The open-class measures ride along with the CLOSED ones: the reference
    // assembles X and Q only when the dispatched method produced them as a
    // by-product (which is `mva` alone) and returns both empty otherwise, so a
    // caller can test emptiness to decide whether it must derive the measures
    // itself. Attaching Qopen to an otherwise empty result would answer that
    // test wrongly and silence the caller's own open-chain formula.
    bool have_measures = false;
    const auto attach_open = [&]() {
        if (ocl.empty() || !have_measures) return;
        res.Q = Matrix<T>(M0, R, zero);
        for (std::size_t i = 0; i < M0; ++i)
            for (std::size_t r : ocl) res.Q(i, r) = Qopen(i, r);
        res.X.assign(R, zero);
        for (std::size_t r : ocl) res.X[r] = lam[r];
    };
    const auto finish = [&](const T& gcore) {
        res.G = Gscale * gcore;
        res.lG = num_traits<T>::log_as_double(res.G);
        attach_open();
        return res;
    };
    // The estimator finish: lG is the answer and G is derived from it, so a
    // constant past the double range still reports a usable logarithm.
    const auto finish_log = [&](double lgcore) {
        res.lG = lGscale + lgcore;
        res.G = num_traits<T>::from_double(std::exp(res.lG));
        attach_open();
        return res;
    };

    // ---- degenerate cases, in closed form ------------------------------------
    T Lsum = zero;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < R1; ++k) Lsum += L2(i, k);

    if (M == 0 || !(Lsum > atol)) {
        // All demands zero: the whole population sits in the delay.
        std::vector<std::size_t> all(R1);
        for (std::size_t k = 0; k < R1; ++k) all[k] = k;
        return finish(Ztot > atol ? delayG(all) : one);
    }
    if (M == 1 && !(Ztot > atol)) {
        // Single station, no delay: G = (sum N)! / prod N_r! * prod L_r^{N_r}.
        long tot = 0;
        for (int v : N2) tot += v;
        T g = num_factorial<T>(static_cast<unsigned>(tot));
        for (std::size_t k = 0; k < R1; ++k)
            g *= num_pow_int(L2(0, k), static_cast<unsigned>(N2[k])) /
                 num_factorial<T>(static_cast<unsigned>(N2[k]));
        return finish(g);
    }
    if (!(Ztot > atol)) {
        // M identical replicas, no delay: the multiset-count closed form.
        bool identical = true;
        for (std::size_t i = 1; i < M && identical; ++i)
            for (std::size_t k = 0; k < R1; ++k)
                if (L2(i, k) != L2(0, k)) {
                    identical = false;
                    break;
                }
        if (identical) {
            long tot = 0;
            for (int v : N2) tot += v;
            T g = num_factorial<T>(static_cast<unsigned>(tot + M - 1)) /
                  num_factorial<T>(static_cast<unsigned>(M - 1));
            for (std::size_t k = 0; k < R1; ++k)
                g *= num_pow_int(L2(0, k), static_cast<unsigned>(N2[k])) /
                     num_factorial<T>(static_cast<unsigned>(N2[k]));
            return finish(g);
        }
    }

    // ---- classes whose jobs never leave the delay -----------------------------
    std::vector<std::size_t> zdem, nzdem;
    for (std::size_t k = 0; k < R1; ++k) {
        T s = zero;
        for (std::size_t i = 0; i < M; ++i) s += L2(i, k);
        (s > atol ? nzdem : zdem).push_back(k);
    }
    const T Gzdem = zdem.empty() ? one : delayG(zdem);

    Matrix<T> L3(M, nzdem.size());
    Matrix<T> Z3(Zag.rows(), nzdem.size());
    std::vector<int> N3(nzdem.size(), 0);
    for (std::size_t a = 0; a < nzdem.size(); ++a) {
        for (std::size_t i = 0; i < M; ++i) L3(i, a) = L2(i, nzdem[a]);
        for (std::size_t i = 0; i < Zag.rows(); ++i) Z3(i, a) = Zag(i, nzdem[a]);
        N3[a] = N2[nzdem[a]];
    }

    // Mean values produced as a by-product come back in the REDUCED, SCALED and
    // REORDERED problem: demand-bearing stations only, nonzero-demand classes
    // only, on demands divided by scalevec. Undo all three, or the caller reads
    // a permuted throughput of the wrong magnitude.
    const auto attach_measures = [&](const std::vector<T>& Xr, const Matrix<T>& Qr) {
        res.X.assign(R, num_traits<T>::from_int(0));
        res.Q = Matrix<T>(M0, R, num_traits<T>::from_int(0));
        for (std::size_t a = 0; a < nzdem.size() && a < Xr.size(); ++a) {
            const std::size_t r = nnz[nzdem[a]];
            res.X[r] = Xr[a] / scalevec[nzdem[a]];
            for (std::size_t i = 0; i < M && i < Qr.rows(); ++i) res.Q(demSt[i], r) = Qr(i, a);
        }
    };

    const std::size_t Rc = nzdem.size();
    T Z3tot = zero;
    for (std::size_t i = 0; i < Z3.rows(); ++i)
        for (std::size_t a = 0; a < Rc; ++a) Z3tot += Z3(i, a);

    // Chain-level views the estimators want: Z summed over its rows, and the
    // population as field elements rather than ints.
    std::vector<T> Zv(Rc, zero);
    for (std::size_t i = 0; i < Z3.rows(); ++i)
        for (std::size_t a = 0; a < Rc; ++a) Zv[a] += Z3(i, a);
    std::vector<T> Nv(Rc, zero);
    for (std::size_t a = 0; a < Rc; ++a) Nv[a] = num_traits<T>::from_int(N3[a]);
    long Nsum3 = 0;
    for (int v : N3) Nsum3 += v;

    // ---- the estimator ladder -------------------------------------------------
    // These are the branches of `compute_norm_const` that answer with a log:
    // asymptotic expansions, quadratures and Monte Carlo estimators. They are
    // separated from the exact switch below because they must return through
    // finish_log, and because the whole group is meaningless -- not merely
    // inaccurate -- in an exact field, where it is refused by name.
    const bool default_multi =
        (method == NcMethod::Default || method == NcMethod::Adaptive) && M > 1;
    const bool default_big_repairman = (method == NcMethod::Default ||
                                        method == NcMethod::Adaptive) &&
                                       M == 1 && Z3tot > atol && Nsum3 >= 10000;
    const bool estimator =
        default_multi || default_big_repairman || method == NcMethod::Clw ||
        method == NcMethod::Cub || method == NcMethod::Gm || method == NcMethod::Kt ||
        method == NcMethod::Bkt || method == NcMethod::Lekt ||
        method == NcMethod::Bk || method == NcMethod::Bkue ||
        method == NcMethod::Lc || method == NcMethod::LcUe ||
        method == NcMethod::Le || method == NcMethod::Ble || method == NcMethod::Ls ||
        method == NcMethod::Aghq ||
        method == NcMethod::Is ||
        method == NcMethod::Mci || method == NcMethod::Imci ||
        method == NcMethod::Mcmc || method == NcMethod::Sampling ||
        method == NcMethod::Mmint2 || method == NcMethod::Gleint ||
        method == NcMethod::Pana || method == NcMethod::Propfair ||
        // 'rgf' answers with a log in both arities: the grouped convolution at
        // one class, the residue recursion of Harrison-Coury Thm 1 beyond it.
        method == NcMethod::Rgf ||
        // 'divdiff' evaluates alternating sums as signed log-sum-exps, so it
        // answers with a log like the estimators do, exact though it is.
        method == NcMethod::Divdiff;
    if (estimator) {
        if constexpr (!num_traits<T>::has_transcendental) {
            pfqn_nc_refuse(default_multi || default_big_repairman
                               ? std::string("default (the multi-station cub / le branch)")
                               : std::string(nc_method_name(method)));
        } else {
            const double lgz = num_traits<T>::log_as_double(Gzdem);
            McRng rng(static_cast<std::uint64_t>(nopt.seed));
            const T fineTol = num_traits<T>::from_double(1e-8);  // GlobalConstants.FineTol
            // The order the reference raises as far as a fixed cost budget allows.
            const auto cub_budget_order = [&]() {
                const double Cmax = static_cast<double>(M * Rc) * 125000.0;
                const int maxorder =
                    static_cast<int>(std::min<double>(std::ceil((Nsum3 - 1) / 2.0), 16.0));
                double tot = 0.0;
                int order = 0;
                while (order < maxorder) {
                    const double next =
                        static_cast<double>(Rc) *
                        nck(static_cast<int>(M) + 2 * (order + 1), static_cast<int>(M) - 1);
                    if (tot + next > Cmax) break;
                    ++order;
                    tot += next;
                }
                // Cmax prices neither the Grundmann-Moeller node count nor the
                // think-time v-integration, so the order is re-priced against
                // the true evaluation count and lowered until it fits.
                double Zsum = 0.0;
                for (std::size_t a = 0; a < Rc; ++a) Zsum += num_traits<T>::to_double(Zv[a]);
                while (order > 0 &&
                       pfqn_cub_evals(static_cast<int>(M), order, Zsum) > CUB_MAX_EVALS)
                    --order;
                return order;
            };
            switch (method) {
                case NcMethod::Default:
                case NcMethod::Adaptive:
                    // ONE ESTIMATOR ANSWERS THE WHOLE FAMILY. The
                    // divided-difference closed form of Casale (SIGMETRICS 2017)
                    // is exact here and was briefly tried first on
                    // M>1 && Rc==1 && sum(Z)==0, but the default route does not
                    // serve a single constant: the analyzer differences it at
                    // N-e_r for X and at the AUGMENTED shape for Q, one extra
                    // class holding one job at station i. That shape has Rc+1
                    // classes, which the closed form refuses at any sizeable
                    // population (the outer sum's cancellation), so it kept the
                    // cubature while G(N) turned exact. Mixing the two costs more
                    // than either: on mqn_singleserver_ps the closed-form G(N)
                    // under cubature numerators left sum_i Q_i at 99.500 of
                    // N=100, and the conservation rescale then moved the entire
                    // cubature error into X, 0.5% against the 0.06% the cubature
                    // ratio carries on its own. 'divdiff' stays a NAMED method,
                    // where the caller owns the whole family.
                    if (M > 1 && Nsum3 < 1000) {
                        res.method = "cub";
                        return finish_log(
                            lgz + num_traits<T>::to_double(
                                      pfqn_cub(L3, N3, Zv, cub_budget_order(), fineTol).lG));
                    }
                    // BLE on the default path: strictly better on lG and it
                    // cancels in G(N-e_r)/G(N). "le" stays the published form.
                    res.method = "ble";
                    // Birman-Kogan Algorithm 2 supplies the MEAN VALUES here.
                    // The caller's fallback differences lG at R+M*R reduced
                    // populations, which on many stations is both dearer and
                    // ~300x less accurate than the load concealment fixed point. Gated
                    // on the station count, since load concealment is mean
                    // field in M: see _kb/06-solver-catalog.md.
                    if (L3.rows() >= 10 && L3.cols() > 1) {
                        bool closed = !Nv.empty();
                        T tot = num_traits<T>::from_int(0);
                        for (std::size_t r = 0; r < Nv.size(); ++r) {
                            if (Nv[r] < num_traits<T>::from_int(0)) closed = false;
                            tot = tot + Nv[r];
                        }
                        if (closed && tot > num_traits<T>::from_int(0)) {
                            BkLcResult<T> thin = pfqn_bklc(L3, Nv, Zv, "mva", 1e-10, 1000);
                            attach_measures(thin.X, thin.Q);
                            res.method = "ble/lc";
                        }
                    }
                    return finish_log(lgz +
                                      num_traits<T>::to_double(pfqn_ble(L3, Nv, Zv).lG));
                case NcMethod::Clw:
                    res.method = "clw";
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_clw(L3, N3, Zv).lG));
                case NcMethod::Cub:
                case NcMethod::Gm: {
                    // The reference's exact order for the Z = 0 branch.
                    const int order = static_cast<int>(std::ceil((Nsum3 - 1) / 2.0));
                    res.method = nc_method_name(method);
                    return finish_log(
                        lgz + num_traits<T>::to_double(pfqn_cub(L3, N3, Zv, order, fineTol).lG));
                }
                case NcMethod::Kt:
                    res.method = "kt";
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_kt(L3, Nv, Zv).lG));
                case NcMethod::Bkt:
                    // KT minus the exact Stirling remainder of each Laplaced class
                    res.method = "bkt";
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_bkt(L3, Nv, Zv).lG));
                case NcMethod::Lekt:
                    // the estimator ble and bkt both compute, on the cheaper side
                    res.method = "lekt";
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_lekt(L3, Nv, Zv).lG));
                case NcMethod::Bk:
                    res.method = "bk";
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_bk(L3, Nv, Zv).lG));
                case NcMethod::Bkue: {
                    // The uniform expansion is single chain by construction; the
                    // multichain fallback is the saddle point of the same paper,
                    // which is also how the analyzer reaches this branch, since it
                    // conditions on a station population with an auxiliary class.
                    if (L3.cols() > 1) {
                        res.method = "bkue/bk";
                        return finish_log(lgz + num_traits<T>::to_double(pfqn_bk(L3, Nv, Zv).lG));
                    }
                    std::vector<T> Dv(L3.rows());
                    for (std::size_t i = 0; i < L3.rows(); ++i) Dv[i] = L3(i, 0);
                    res.method = "bkue";
                    return finish_log(
                        lgz + num_traits<T>::to_double(
                                  pfqn_bkue(Dv, Nv[0], Zv.empty() ? num_traits<T>::from_int(0) : Zv[0]).lG));
                }
                case NcMethod::Lc:
                case NcMethod::LcUe: {
                    // Algorithm 2 returns mean values, not a multichain constant;
                    // the saddle point that seeds it supplies lG on the same
                    // asymptotics. The fixed point converges linearly and slowly,
                    // so it is iterated to the method's own accuracy rather than to
                    // a solver-level reporting tolerance.
                    BkLcResult<T> thin = pfqn_bklc(
                        L3, Nv, Zv, method == NcMethod::LcUe ? "ue" : "mva", 1e-10, 1000);
                    attach_measures(thin.X, thin.Q);
                    res.method = nc_method_name(method);
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_bk(L3, Nv, Zv).lG));
                }
                case NcMethod::Le:
                    res.method = "le";
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_le(L3, Nv, Zv).lG));
                case NcMethod::Ble:
                    // LE plus the empirical eps->0 correction; see _kb/03-api-layer.md
                    res.method = "ble";
                    return finish_log(lgz + num_traits<T>::to_double(pfqn_ble(L3, Nv, Zv).lG));
                case NcMethod::Aghq:
                    // adaptive Gauss-Hermite over the simplex; q=1 would be "le";
                    // options.config.aghq_nodes overrides the node count.
                    res.method = "aghq";
                    return finish_log(lgz + num_traits<T>::to_double(
                        pfqn_aghq(L3, Nv, Zv,
                                  nopt.aghq_nodes < 1 ? 1 : nopt.aghq_nodes).lG));
                case NcMethod::Ls:
                    res.method = "ls";
                    return finish_log(
                        lgz + num_traits<T>::to_double(pfqn_ls(L3, Nv, Zv, nopt.samples, rng).lG));
                case NcMethod::Is:
                    res.method = "is";
                    return finish_log(lgz + pfqn_is(L3, N3, Zv, nopt.samples, rng).lG);
                case NcMethod::Mci:
                    res.method = "mci";
                    return finish_log(
                        lgz + pfqn_mci(L3, N3, Zv, nopt.samples, MciVariant::Mci, rng).lG);
                case NcMethod::Imci:
                    res.method = "imci";
                    return finish_log(
                        lgz + pfqn_mci(L3, N3, Zv, nopt.samples, MciVariant::Imci, rng).lG);
                case NcMethod::Mcmc: {
                    // Chen-O'Cinneide REGULARIZATION (TOMACS 8(3), 1998). The chain is
                    // simulated on the regularized network, which shares the steady-state
                    // distribution of the original one, so X and Q come back through the
                    // X/Q channel like 'lc'. What it does NOT return is the constant
                    // itself: the algorithm estimates the RATIOS G(N-e_r)/G(N), never G,
                    // so the lG below is the BLE expansion and is not part of the paper.
                    // It cancels out of every mean value the analyzer reports; only
                    // getProbNormConstAggr reads it.
                    const McmcResult<T> mc = pfqn_mcmc(L3, N3, Zv, std::vector<double>(),
                                                       nopt.samples, nopt.mcmc_batches,
                                                       nopt.mcmc_burnin, rng);
                    attach_measures(mc.X, mc.Q);
                    res.method = "mcmc";
                    if (L3.rows() > 1)
                        return finish_log(lgz + num_traits<T>::to_double(pfqn_ble(L3, Nv, Zv).lG));
                    return finish_log(lgz + pfqn_comomrm(L3, N3, Z3, 1).lG);
                }
                case NcMethod::Sampling:
                    // The reference picks by shape: one station is a repairman
                    // integral, more stations than classes favour the Monte
                    // Carlo integral, and otherwise the logistic sampler.
                    if (M == 1) {
                        res.method = "sampling";
                        return finish_log(
                            lgz + pfqn_mmsample2(L3, N3, Zv, nopt.samples, rng).lG);
                    }
                    if (M > Rc) {
                        res.method = "imci";
                        return finish_log(
                            lgz + pfqn_mci(L3, N3, Zv, nopt.samples, MciVariant::Imci, rng).lG);
                    }
                    res.method = "ls";
                    return finish_log(
                        lgz + num_traits<T>::to_double(pfqn_ls(L3, Nv, Zv, nopt.samples, rng).lG));
                case NcMethod::Mmint2:
                case NcMethod::Gleint: {
                    if (M > 1) {
                        // The reference warns and returns lG = []; see `valid`.
                        res.method = nc_method_name(method);
                        res.valid = false;
                        res.lG = 0.0;
                        res.G = zero;
                        return res;
                    }
                    std::vector<T> Lrow(Rc, zero);
                    for (std::size_t a = 0; a < Rc; ++a) Lrow[a] = L3(0, a);
                    res.method = nc_method_name(method);
                    return finish_log(
                        lgz + num_traits<T>::to_double(
                                  pfqn_mmint2_gausslegendre(Lrow, Nv, Zv).lG));
                }
                case NcMethod::Pana: {
                    const PanaceaResult<T> pa = pfqn_panacea(L3, N3, Zv);
                    res.method = "pana";
                    if (!pa.normalUsage) {
                        // The reference warns and returns lG = []; see the
                        // `valid` field. Deliberate, user-ruled 2026-07-25.
                        res.valid = false;
                        res.lG = 0.0;
                        res.G = zero;
                        return res;
                    }
                    return finish_log(lgz + num_traits<T>::to_double(pa.lG));
                }
                case NcMethod::Propfair:
                    res.method = "propfair";
                    return finish_log(lgz +
                                      num_traits<T>::to_double(pfqn_propfair(L3, Nv, Zv).lG));
                case NcMethod::Divdiff: {
                    // Divided-difference closed form, Eqs. (15) and (16).
                    // Load-independent single-server queues only: a think time
                    // needs the integral form of Corollary 3.4, which is not
                    // implemented. Unlike the default route this one keeps
                    // whatever the expression returns, since a caller that named
                    // the method has no fallback.
                    if (Z3tot > zero)
                        throw InputError(
                            "pfqn_nc: the 'divdiff' method requires a model without think time, "
                            "which needs the integral form of Corollary 3.4. Use 'ca' or "
                            "'default'");
                    const ExplicitResult<T> ex = pfqn_explicit(L3, N3);
                    res.method = "divdiff/" + ex.method;
                    return finish_log(lgz + num_traits<T>::to_double(ex.lG));
                }
                case NcMethod::Rgf: {
                    if (Rc == 1) {
                        std::vector<T> Lcol(M, zero);
                        for (std::size_t i = 0; i < M; ++i) Lcol[i] = L3(i, 0);
                        res.method = "rgf";
                        return finish_log(
                            lgz + num_traits<T>::to_double(pfqn_rgf(Lcol, N3[0], Zv[0]).lG));
                    }
                    // Multiclass: the residue recursion, with think times carried
                    // by the Bertozzi-McKenna truncation neither RGF paper has.
                    // That sum is ALTERNATING, so pfqn_rgfmc refuses rather than
                    // return a wrong lG; the exact convolution answers those and
                    // the reported method says so.
                    try {
                        const T lg = pfqn_rgfmc<T>(L3, N3, Zv).lG;
                        res.method = "rgf";
                        return finish_log(lgz + num_traits<T>::to_double(lg));
                    } catch (const InputError&) {
                        Matrix<T> Zm(1, Rc);
                        for (std::size_t r = 0; r < Rc; ++r) Zm(0, r) = Zv[r];
                        res.method = "rgf/ca";
                        return finish_log(
                            lgz + num_traits<T>::to_double(pfqn_ca<T>(L3, N3, Zm).lG));
                    }
                }
                default:
                    break;
            }
        }
    }

    // ---- dispatch -------------------------------------------------------------
    T gcore = one;
    switch (method) {
        case NcMethod::Ca:
            gcore = pfqn_ca(L3, N3, Z3).G;
            res.method = "ca";
            break;
        case NcMethod::Rgf:
            // Multiclass only: the single-class recursion returned above. The
            // aggregate marginals reach this branch too, since the analyzer
            // augments a single-class model with an auxiliary class to
            // condition on a station population.
            gcore = pfqn_ca(L3, N3, Z3).G;
            res.method = "rgf/ca";
            break;
        case NcMethod::Ger:
            // Residue closed form of the same generating function 'clw' inverts
            // numerically. It sits in the EXACT switch, not among the estimators,
            // because every operation is a product, a quotient or a binomial
            // coefficient: in the exact backend it agrees with pfqn_ca bit for bit.
            // A class eliminated by residues enters only as a pole ORDER, so its
            // population is free; the term count grows as C(S+M-1,M-1) per further
            // elimination instead, and the maxterms cap REFUSES rather than
            // truncating.
            gcore = pfqn_gerasimov(L3, N3, Z3).G;
            res.method = "ger";
            break;
        case NcMethod::Mva: {
            // MVA by-product means rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
            const MvaResult<T> mva = pfqn_mva(L3, N3, Z3);
            gcore = mva.G;
            res.X.assign(R, zero);
            res.Q = Matrix<T>(M0, R, zero);
            for (std::size_t a = 0; a < Rc; ++a) {
                const std::size_t r = nnz[nzdem[a]];
                res.X[r] = mva.XN[a] / scalevec[nzdem[a]];
                for (std::size_t i = 0; i < M; ++i) res.Q(demSt[i], r) = mva.QN(i, a);
            }
            res.method = "mva";
            have_measures = true;
            break;
        }
        case NcMethod::Recal:
            if (Z3tot > atol)
                throw UnsupportedError(
                    "pfqn_nc: RECAL is available only for models with zero think time; this model "
                    "has a delay");
            gcore = pfqn_recal(L3, N3, Z3).G;
            res.method = "recal";
            break;
        case NcMethod::Exact: {
            // RECAL routing rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
            long tot = 0;
            for (int v : N3) tot += v;
            if (M >= Rc || tot > 10 || Z3tot > atol) {
                gcore = pfqn_ca(L3, N3, Z3).G;
                res.method = "exact/ca";
            } else {
                gcore = pfqn_recal(L3, N3, Z3).G;
                res.method = "exact/recal";
            }
            break;
        }
        case NcMethod::Comom:
            // THIS THROW IS THE REFERENCE'S OWN, and it sits a few lines from
            // branches that DECLINE by returning an empty lG (see `valid`).
            // The rationale is quoted verbatim from pfqn_nc.m:278-281, so that
            // a reader who sees zeros next door does not "harmonise" it away:
            //
            //   "This used to emit a warning gated on options.verbose and then
            //    return lG = [], which the caller turned into all-zero queue
            //    lengths while still reporting a completed analysis. A silently
            //    zeroed result is worse than no result: refuse instead."
            //
            // Declining and refusing are different, and the reference does both
            // on purpose. The 2026-07-25 empty-result ruling covers the first.
            if (Rc > 1 && M > 1)
                throw InputError(
                    "pfqn_nc: the 'comom' method supports a single queueing station, but this "
                    "model has more. Use 'default', 'ca' or 'exact'.");
            if (Rc > 1) {
                gcore = pfqn_comomrm(L3, N3, Z3).G;
                res.method = "comom";
            } else {
                gcore = pfqn_ca(L3, N3, Z3).G;
                res.method = "ca";
            }
            break;
        case NcMethod::Default:
        case NcMethod::Adaptive:
            if (M > 1) {
                // unported cub/le routing rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                pfqn_nc_refuse("default (multi-station branch: cub / le)");
            } else if (!(Z3tot > atol)) {
                // REFERENCE DEFECT (single-queue-no-delay lG): see _kb/03-api-layer.md (cpp port notes: pfqn)
                long tot = 0;
                for (int v : N3) tot += v;
                T g = num_factorial<T>(static_cast<unsigned>(tot));
                for (std::size_t a = 0; a < Rc; ++a)
                    g *= num_pow_int(L3(0, a), static_cast<unsigned>(N3[a])) /
                         num_factorial<T>(static_cast<unsigned>(N3[a]));
                gcore = g;
                res.method = "exact";
            } else {
                gcore = pfqn_comomrm(L3, N3, Z3).G;
                res.method = "comom";
            }
            break;
        default:
            // Every remaining name is an estimator and returned above.
            pfqn_nc_refuse(nc_method_name(method));
    }

    return finish(Gzdem * gcore);
}

/** Overload with the reference's default sample count, seed and tolerance. */
template <class T>
NcDispatchResult<T> pfqn_nc(const std::vector<T>& lambda, const Matrix<T>& L,
                            const std::vector<int>& N, const Matrix<T>& Z, NcMethod method,
                            const T& atol) {
    return pfqn_nc(lambda, L, N, Z, method, atol, NcOptions());
}

/** Overload with the exact (zero-tolerance) filters and no open classes. */
template <class T>
NcDispatchResult<T> pfqn_nc(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                            NcMethod method) {
    return pfqn_nc(std::vector<T>(), L, N, Z, method, num_traits<T>::from_int(0), NcOptions());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_NC_H
