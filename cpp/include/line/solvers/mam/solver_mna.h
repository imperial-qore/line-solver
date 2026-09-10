/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MNA_H
#define LINE_SOLVERS_MAM_SOLVER_MNA_H

/**
 * Port of `solver_mna_open.m` and `solver_mna_closed.m`, the two analyzers
 * behind SolverMAM's `mna` method.
 *
 * THE METHOD. MNA is QNA's flow decomposition with QNA's isolated-station
 * solution REPLACED, at every FCFS station, by a matrix-analytic one. Each
 * sweep superposes the per-class flows into a station as a rate a1(i,r) and a
 * squared coefficient of variation a2(i,r), solves the station in isolation,
 * and splits the departure stream along the outgoing arcs with the exact
 * Bernoulli-thinning rule f2 = 1 + p (d2 - 1). What makes it MAM rather than
 * QNA is the last step: the converged (a1, a2) pair is turned back into a
 * phase-type arrival process -- one APH per class fitted to those two moments,
 * marked and superposed into an MMAP -- and the station's queue length comes
 * from MMAP[K]/PH[K]/1 FCFS rather than from a Whitt waiting-time formula. The
 * flow SCVs d2 that drive the NEXT sweep still come from QNA's departure
 * formula, so the fixed point is QNA's and only the reported queue lengths are
 * matrix-analytic.
 *
 * THE TWO ANALYZERS DIFFER IN WHAT DRIVES THE OUTER LOOP, not in the sweep.
 * The open one has its arrival rates fixed by the sources and runs a single
 * flow fixed point. The closed one has no source: it wraps the same flow fixed
 * point in an OUTER BISECTION on the per-class throughput, bracketed below by 0
 * and above by the slowest service rate over the finite-server stations, and
 * driven against the population target N. That is why the closed analyzer
 * evaluates the matrix-analytic station solve once per outer step, and why it
 * asks for the queue-length DISTRIBUTION (truncated at the class population)
 * where the open one asks only for the mean.
 *
 * WHAT THE REFERENCE GETS WRONG, REPRODUCED RATHER THAN CORRECTED. Three
 * things, each of which changes reported numbers:
 *
 *  - A PS STATION IN THE OPEN ANALYZER REPORTS NOTHING. `solver_mna_open.m`'s
 *    PS branch assigns to `TN`, `UN`, `QN`, `RN`, which are fresh undefined
 *    variables in that function -- the metrics it returns are `T`, `U`, `Q`,
 *    `R`. So a PS station keeps the zeros it was initialised with, and its
 *    departure SCV d2 is never set either. The closed analyzer's PS branch
 *    writes the right names and does work. Reproduced: silently redirecting
 *    the writes would report queue lengths MATLAB does not report.
 *  - THE CLOSED ANALYZER READS ONE CLASS'S DISTRIBUTION FOR EVERY CLASS.
 *    `[pdistr] = MMAPPH1FCFS(..., 'ncDistr', maxLevel)` captures only the FIRST
 *    output, i.e. class 1's marginal, and the per-class truncation loop then
 *    truncates THAT at each class's own population. The reference labels the
 *    block "rough approximation" itself.
 *  - THROUGHPUT IS NEVER REPORTED. Both analyzers initialise `X = zeros(1,K)`
 *    and never assign it, so `getAvg`'s per-class throughput column is zero
 *    however the model is solved. Station throughputs in `T` are correct; it is
 *    only the class-level `X` that is dead.
 *
 * ONE REFERENCE BRANCH IS REFUSED RATHER THAN REPRODUCED. `config.dep_scv =
 * 'etaqa'` is dead in MATLAB -- `qbd_depproc_jointmom` raises on every input it
 * can be handed, so the try/catch around it always takes the QNA fallback --
 * while the ported `qbd_depproc_jointmom` works. Answering with the working one
 * would report a departure SCV the reference never produces, so the option is
 * refused by name instead; see `solver_mna_open`.
 *
 * REFERENCE INDEXING, CHECKED RATHER THAN ASSUMED, exactly as `solver_qna.h`
 * does for the same reason: both files index the STATEFUL-indexed `sn.rt` with
 * station indices, and the closed one additionally indexes the CLASS-indexed
 * `sn.njobs` and the length-C `lambda` with the same running index. Both are
 * silently correct only on the model shapes the method is advertised for, and
 * are refused by name otherwise.
 *
 * SELF-LOOPING CLASSES. `sn.isslc` guards three blocks in the closed analyzer.
 * The C++ `JobClassType` is OPEN or CLOSED only, so no model this port can
 * build enters them and they are not transcribed.
 *
 * Arithmetic: DOUBLE-ONLY, gated the way `solver_mam_basic.h` gates it. The
 * flow fixed point stops on a tolerance, the APH fit takes a ceiling of a real
 * reciprocal, and MMAP[K]/PH[K]/1 runs the ADDA doubling iteration.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "line/api/da/da_fpi.h"
#include "line/api/da/da_traffic_superpos.h"
#include "line/api/mam/aph_fit_moments.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_assemble.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/api/npfqn/npfqn_traffic_split_rr.h"
#include "line/api/qsys/qsys_mmck.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/solvers/mam/solver_mam_bmap.h"  // mam_detect_mmck, shared with solver_mam_basic
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * The `options.config` fields the two MNA analyzers read.
 *
 * Kept out of `MamOptions` because neither field is a SolverMAM option: both
 * are written by the analyzer itself onto a LOCAL copy of `options.config`, so
 * a caller has no way to reach them through the solver's option surface.
 */
struct MnaConfig {
    /**
     * `config.dep_scv`, read only by the OPEN analyzer. 'qna' takes Whitt's
     * departure-SCV formula. 'etaqa' is refused by name; see the branch.
     */
    std::string dep_scv = "qna";
};

namespace mna_detail {

using lang::GlobalConstants;
using lang::ProcessType;
using lang::SchedStrategy;

// aph_from_2moments and aph_fit_mean_scv now live in
// api/mam/aph_fit_moments.h: the closed setup/delay-off branch of
// solver_mam_basic needs the same APH.fitMeanAndSCV, and that header
// already includes this one, so keeping them here would have been a
// cycle. Re-exported so every existing caller is unchanged.
using mam::aph_from_2moments;
using mam::aph_fit_mean_scv;

// `mam_detect_mmck` used to be transcribed a second time here, under the name
// `detect_mmck` and a 0-based bool/out-param signature, with a comment saying
// "the two must stay in step". Collapsed onto the single definition in
// solver_mam_bmap.h: a duplicate whose own comment admits it must be kept in
// step is the defect class that already cost this tree once, in the two live
// ports of solver_fluid_initsol.m. The surviving form is the 1-based
// struct-returning one because it is the one carrying a station-index range
// check; the bodies were otherwise character-for-character the same.

/** The gates both analyzers share, named after the indexing they protect. */
template <class T>
void check_gates(const qn::NetworkStruct<T>& L, const std::string& who) {
    // The Fork test comes first because a Fork-Join model also fails the
    // stateful-count test, and the reference's own refusal is the Fork one.
    if (L.has_fork())
        throw UnsupportedError(who +
                               ": Fork nodes are not supported yet by the QNA-family solvers, as "
                               "the reference's line_error states");
    if (L.nof_stateful() != L.nstations)
        throw UnsupportedError(
            who + ": the reference indexes the stateful-indexed sn.rt with station indices, which "
                  "is only correct when every stateful node is a station; this model has " +
            std::to_string(L.nof_stateful()) + " stateful nodes and " +
            std::to_string(L.nstations) + " stations");
}

/** `sn.pie` and `sn.proc{ist}{k}{1}` at every station the analyzers solve. */
template <class T>
std::vector<std::vector<PhService<T>>> service_laws(const qn::NetworkStruct<T>& L) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    std::vector<std::vector<PhService<T>>> svc(M, std::vector<PhService<T>>(K));
    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy sc = L.stations[i].sched;
        if (!(sc == SchedStrategy::FCFS || sc == SchedStrategy::INF || sc == SchedStrategy::PS))
            continue;
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.has_service_law(i, r) || !(L.rates(i, r) > zero)) {
                // The reference's `any(isnan(D0))` guard: a class this station
                // never serves gets an Immediate service rather than a NaN pair.
                // has_service_law also covers the Join, whose rates are Inf but
                // whose process is the NaN Coxian this guard exists for.
                const T imm = num_traits<T>::from_double(GlobalConstants::Immediate);
                svc[i][r].sigma.assign(1, one);
                svc[i][r].S = Matrix<T>(1, 1, T(-imm));
                continue;
            }
            const Map<T> ph = lang::dist_to_map(L.service[i][r]);
            svc[i][r].sigma = map_pie(ph);
            svc[i][r].S = ph.D0;
        }
    }
    return svc;
}

// `station_visits` now lives in basic_detail (solver_mam_basic.h), for the same
// reason `mam_detect_mmck` does: solver_mam_basic and solver_mam_basic_mmap need
// the identical cellsum(sn.visits), and a third body to keep in step is the
// defect class this file already collapsed once. Re-exported so every caller
// below is unchanged.
using basic_detail::station_visits;

/**
 * f2 on every arc that does not end at a Source, the shared initial state.
 *
 * `kRR` is the round-robin split degree of each station-class
 * (`npfqn_traffic_split_rr`); it is all ones in the closed analyzer, which the
 * reference does not correct, and the entry then reduces to the plain 1 the
 * Bernoulli-thinning rule gives at d2 = 1.
 */
template <class T>
Matrix<T> init_flow_scv(const qn::NetworkStruct<T>& L, const std::vector<bool>& is_source,
                        const Matrix<T>& kRR) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    Matrix<T> f2(M * K, M * K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            if (is_source[j]) continue;
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t s = 0; s < K; ++s) {
                    const T p = L.rt(i * K + r, j * K + s);
                    if (p > zero) f2(i * K + r, j * K + s) = T(one + p * T(one - kRR(i, r)));
                }
        }
    return f2;
}

/**
 * The superposition step, identical in both analyzers: a1 accumulates the
 * routed throughput and a2 the rate-weighted mixture of the incoming flow SCVs.
 */
template <class T>
void superpose(const qn::NetworkStruct<T>& L, const Matrix<T>& Tp, const Matrix<T>& f2,
               Matrix<T>& a1, Matrix<T>& a2) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    for (std::size_t i = 0; i < M; ++i) {
        T lambda_i = zero;
        for (std::size_t k = 0; k < K; ++k) lambda_i += Tp(i, k);
        for (std::size_t r = 0; r < K; ++r) {
            a1(i, r) = zero;
            a2(i, r) = zero;
        }
        for (std::size_t j = 0; j < M; ++j)
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t s = 0; s < K; ++s) {
                    const T p = L.rt(j * K + s, i * K + r);
                    if (!(p > zero)) continue;
                    a1(i, r) = T(a1(i, r) + Tp(j, s) * p);
                    // A station with no throughput divides by zero in MATLAB and
                    // carries the resulting NaN into a2; the guard keeps a2 at
                    // zero there, which is the value every downstream branch
                    // reads once the NaN sweep at the end has run.
                    if (lambda_i > zero)
                        a2(i, r) = T(a2(i, r) + T(one / lambda_i) * f2(j * K + s, i * K + r) *
                                                    Tp(j, s) * p);
                }
    }
}

/**
 * The splitting step, identical in both analyzers.
 *
 * A flow carrying a fraction p of a stream dispatched one-in-k is the k-fold
 * convolution thinned at q = k p, hence C^2 = 1 + p (d2 - k); k = 1 is the
 * Bernoulli-thinned renewal stream and the only case the closed analyzer sees.
 */
template <class T>
void split(const qn::NetworkStruct<T>& L, const std::vector<bool>& is_source,
           const std::vector<T>& d2, const Matrix<T>& kRR, Matrix<T>& f2) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            if (is_source[j]) continue;
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t s = 0; s < K; ++s) {
                    const T p = L.rt(i * K + r, j * K + s);
                    if (p > zero) f2(i * K + r, j * K + s) = T(one + p * T(d2[i] - kRR(i, r)));
                }
        }
}

/** Whitt's departure-SCV formula, the `dep_scv = 'qna'` branch of both files. */
template <class T>
T qna_departure_scv(const qn::NetworkStruct<T>& L, std::size_t i0, const Matrix<T>& a1,
                    const Matrix<T>& a2, const Matrix<T>& scv, const T& lambda_ist, const T& rho,
                    const T& mi) {
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T mubar = (rho > zero) ? T(lambda_ist / rho) : zero;
    T c2 = T(-one);
    for (std::size_t r = 0; r < L.nclasses; ++r) {
        if (L.disabled[i0][r] || !(L.rates(i0, r) > zero) || !(lambda_ist > zero)) continue;
        const T q = T(mubar / mi / L.rates(i0, r));
        c2 += T(T(a1(i0, r) / lambda_ist) * q * q * T(scv(i0, r) + one));
    }
    T a2sum = zero;
    for (std::size_t r = 0; r < L.nclasses; ++r) a2sum += a2(i0, r);
    return T(one + T(rho * rho * T(c2 - one) / sqrt(mi)) + T(T(one - rho * rho) * T(a2sum - one)));
}

/**
 * The station's aggregate utilization, as both files compute it: a1 over
 * FineTol + rates, summed over classes and divided by the server count.
 *
 * The FineTol in the denominator is the reference's guard against a zero rate,
 * and it is why a station whose classes are all disabled reports rho = 0 rather
 * than a division by zero.
 */
template <class T>
T station_rho(const qn::NetworkStruct<T>& L, std::size_t i0, const Matrix<T>& a1, const T& mi) {
    const T zero = num_traits<T>::from_int(0);
    const T ftol = num_traits<T>::from_double(GlobalConstants::FineTol);
    T rho = zero;
    for (std::size_t r = 0; r < L.nclasses; ++r) {
        if (L.disabled[i0][r]) continue;  // MATLAB's NaN rate, dropped by isnan
        rho += T(a1(i0, r) / T(ftol + L.rates(i0, r)));
    }
    return T(rho / mi);
}

/**
 * The arrival MMAP a station is solved against: one APH per class fitted to
 * (1/a1, a2), marked as its own class and superposed.
 *
 * A class with no inflow contributes the reference's `map_exponential(Inf)`, a
 * zero-rate order-1 stream. It is carried as a ONE-CLASS MMAP so that the
 * superposition still delivers K marks -- the queue solver downstream is asked
 * for K per-class answers, and a mark-free component would silently shorten
 * that list. `solver_mam_basic.h` reads the same idiom the same way.
 *
 * @param bounded true for the closed analyzer, which superposes through
 *                mmap_super_safe under an order budget; the open one calls the
 *                unbounded mmap_super, its own space_max being dead code
 * @param a1 per-station per-class mean interarrival times; row i0 is read
 * @param a2 per-station per-class interarrival SCVs; row i0 is read
 * @param i0 row index of the station whose arrivals are being built
 * @param K number of classes, and of marks in the returned MMAP
 * @param space_max order budget of the bounded superposition
 */
template <class T>
Mmap<T> arrival_mmap(const Matrix<T>& a1, const Matrix<T>& a2, std::size_t i0, std::size_t K,
                     bool bounded, std::size_t space_max) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Mmap<T> node;
    for (std::size_t k = 0; k < K; ++k) {
        Mmap<T> cur;
        if (a1(i0, k) == zero) {
            cur.D0 = Matrix<T>(1, 1, zero);
            cur.D1 = Matrix<T>(1, 1, zero);
            cur.Dc.assign(1, Matrix<T>(1, 1, zero));
        } else {
            const Map<T> ph = aph_fit_mean_scv(T(one / a1(i0, k)), a2(i0, k));
            cur.D0 = ph.D0;
            cur.D1 = ph.D1;
            cur.Dc.assign(1, ph.D1);
        }
        if (k == 0) node = cur;
        else if (bounded) node = mmap_super_safe(std::vector<Mmap<T>>{node, cur}, space_max);
        else node = mmap_super(node, cur);
    }
    return node;
}

// `zero_nans` moved to basic_detail beside `station_visits`, for the same
// reason: solver_mam_basic_mmap ends with the identical sweep.
using basic_detail::zero_nans;

}  // namespace mna_detail

/**
 * Port of `solver_mna_open.m`.
 *
 * @param L   the refreshed struct; open chains only
 * @param opt SolverMAM's options; `tol` drives the saturation test and doubles
 *            as the fixed point's `iter_tol`, which SolverOptions('MAM') leaves
 *            at the same 1e-4
 * @param cfg the `options.config` fields the analyzer reads
 */
template <class T>
mva::MvaSolution<T> solver_mna_open(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                                    const MnaConfig& cfg = MnaConfig()) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)opt;
        (void)cfg;
        throw UnsupportedError(
            "solver_mna_open: the flow fixed point stops on a tolerance, the APH arrival fit takes "
            "a ceiling of a real reciprocal, and the MMAP[K]/PH[K]/1 station solve runs the ADDA "
            "doubling iteration; rerun this model with --arith double or --arith real");
    } else {
    using namespace mna_detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;
    const T tol = num_traits<T>::from_double(opt.tol);

    check_gates(L, "solver_mna_open");
    for (const qn::JobClass& c : L.classes)
        if (std::isfinite(c.population))
            throw UnsupportedError(
                "solver_mna_open: MNA's open analyzer takes its arrival rates from the class "
                "sources; a closed chain has none and belongs to solver_mna_closed");
    if (cfg.dep_scv == "etaqa")
        throw UnsupportedError(
            "solver_mna_open: config.dep_scv = 'etaqa' reads the departure SCV off "
            "qbd_depproc_jointmom, and the MATLAB routine of that name raises a dimension error "
            "on EVERY input this branch can hand it -- it slices the level-0 vector as pi(1,:) "
            "from a QBD_pi that returns one long row, so v0 is numLevels times too long "
            "(measured: arrival/service orders 1/1, 1/2 and 2/2 all throw). The reference's "
            "try/catch therefore always falls back to the QNA formula, while the ported "
            "qbd_depproc_jointmom takes the correct slice and SUCCEEDS; running it here would "
            "report a departure SCV the reference never produces. Use the default 'qna'");
    if (cfg.dep_scv != "qna")
        throw UnsupportedError("solver_mna_open: unknown config.dep_scv '" + cfg.dep_scv +
                               "'; the reference offers 'qna' and 'etaqa'");

    Matrix<T> S(M, K, zero), scv(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.disabled[i][r] && L.rates(i, r) > zero) S(i, r) = T(one / L.rates(i, r));
            const double v = num_traits<T>::to_double(L.scv(i, r));
            scv(i, r) = std::isnan(v) ? zero : L.scv(i, r);
        }
    const Matrix<T> V = station_visits(L);
    const std::vector<std::vector<PhService<T>>> svc = service_laws(L);

    std::vector<bool> is_source(M, false);
    for (std::size_t i = 0; i < M; ++i)
        is_source[i] = (L.stations[i].nodetype == qn::NodeType::Source);
    // deterministic (round-robin) split degrees, 1 where the split is Markovian
    const Matrix<T> kRR = npfqn::npfqn_traffic_split_rr(L);
    Matrix<T> f2 = init_flow_scv(L, is_source, kRR);

    Matrix<T> Q(M, K, zero), U(M, K, zero), R(M, K, zero), Tp(M, K, zero);
    Matrix<T> a1(M, K, zero), a2(M, K, zero);
    std::vector<T> d2(M, zero), lambda(C, zero);

    // ---- the source streams ----------------------------------------------
    std::vector<T> d2c(C, zero);
    std::size_t sourceIdx = 0;
    for (std::size_t c = 0; c < C; ++c) {
        const std::size_t ref = L.classes[L.inchain[c][0] - 1].refstat;
        sourceIdx = ref;
        std::vector<T> lam_in, scv_in;
        for (std::size_t k : L.inchain[c]) {
            lam_in.push_back(L.disabled[ref - 1][k - 1] ? zero : L.rates(ref - 1, k - 1));
            scv_in.push_back(scv(ref - 1, k - 1));
        }
        for (const T& x : lam_in)
            if (std::isfinite(num_traits<T>::to_double(x))) lambda[c] += x;
        d2c[c] = da::da_traffic_superpos(lam_in, scv_in);
        for (std::size_t a = 0; a < L.inchain[c].size(); ++a)
            Tp(ref - 1, L.inchain[c][a] - 1) = lam_in[a];
    }
    // `d2(sourceIdx) = d2c(sourceIdx,:)*lambda'/sum(lambda)` runs ONCE, after
    // the chain loop, so it seeds only the LAST chain's reference station. The
    // row index into the 1 x C vector d2c is that same station, so MATLAB
    // itself only evaluates this when the station is the first one.
    if (sourceIdx != 1)
        throw UnsupportedError(
            "solver_mna_open: the reference seeds the source departure SCV with "
            "d2c(sourceIdx,:), indexing the 1 x nchains vector d2c by the reference STATION " +
            std::to_string(sourceIdx) +
            "; MATLAB evaluates that only when the reference station is station 1");
    {
        T num = zero, den = zero;
        for (std::size_t c = 0; c < C; ++c) {
            num += T(d2c[c] * lambda[c]);
            den += lambda[c];
        }
        if (den > zero) d2[sourceIdx - 1] = T(num / den);
    }

    // ---- the flow fixed point --------------------------------------------
    auto sweep = [&](const std::vector<T>&,
                     std::size_t itnum) -> std::pair<std::vector<T>, std::vector<T>> {
        std::vector<T> xref(2 * M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k) {
                xref[i * K + k] = a1(i, k);
                xref[M * K + i * K + k] = a2(i, k);
            }

        if (itnum == 1)
            for (std::size_t c = 0; c < C; ++c)
                for (std::size_t m = 0; m < M; ++m)
                    for (std::size_t k : L.inchain[c]) Tp(m, k - 1) = T(V(m, k - 1) * lambda[c]);

        superpose(L, Tp, f2, a1, a2);

        // The reference walks the NODES and maps each to its station; no branch
        // reads another station's iterate, so station order is equivalent.
        for (std::size_t i = 0; i < M; ++i) {
            const SchedStrategy sched = L.stations[i].sched;
            const T mi = num_traits<T>::from_double(L.stations[i].nservers);
            if (sched == SchedStrategy::INF) {
                // MATLAB writes d2(ist,s) = a2(ist,s) across classes but every
                // downstream read is the scalar d2(ist), i.e. column 1.
                d2[i] = a2(i, 0);
                for (std::size_t c = 0; c < C; ++c)
                    for (std::size_t k : L.inchain[c]) {
                        const std::size_t r = k - 1;
                        Tp(i, r) = a1(i, r);
                        U(i, r) = T(S(i, r) * Tp(i, r));
                        Q(i, r) = T(Tp(i, r) * S(i, r) * V(i, r));
                        R(i, r) = (Tp(i, r) > zero) ? T(Q(i, r) / Tp(i, r)) : zero;
                    }
            } else if (sched == SchedStrategy::PS) {
                // Deliberately empty: see the file header. The reference's PS
                // branch writes TN/UN/QN/RN, which are not the metrics it
                // returns, so this station contributes nothing and leaves d2 at
                // whatever the previous sweep left there.
            } else if (sched == SchedStrategy::FCFS) {
                T lambda_ist = zero;
                for (std::size_t r = 0; r < K; ++r) lambda_ist += a1(i, r);
                const T rho = station_rho(L, i, a1, mi);
                if (rho < T(one - tol)) {
                    d2[i] = qna_departure_scv(L, i, a1, a2, scv, lambda_ist, rho, mi);
                } else {
                    for (std::size_t r = 0; r < K; ++r)
                        Q(i, r) = num_traits<T>::from_double(L.classes[r].population);
                    d2[i] = one;
                }
                for (std::size_t r = 0; r < K; ++r) {
                    Tp(i, r) = a1(i, r);
                    U(i, r) = T(Tp(i, r) * S(i, r) / mi);
                }
            } else if (sched != SchedStrategy::EXT) {
                throw UnsupportedError(
                    std::string("solver_mna_open: no isolated-station solution for ") +
                    lang::sched_to_text(sched) + " scheduling at station '" +
                    L.stations[i].name + "'");
            }
        }

        split(L, is_source, d2, kRR, f2);

        std::vector<T> xnew(2 * M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k) {
                xnew[i * K + k] = a1(i, k);
                xnew[M * K + i * K + k] = a2(i, k);
            }
        return std::make_pair(xnew, xref);
    };

    da::FpiOptions fo;
    fo.iter_max = static_cast<std::size_t>(opt.iter_max) + 1;  // the legacy loop ran one extra
    fo.iter_tol = opt.tol;
    fo.nanstop = true;
    const da::FpiResult<T> fr = da::da_fpi<T>(sweep, std::vector<T>(2 * M * K, zero), fo);

    // ---- the matrix-analytic pass over the FCFS stations -------------------
    for (std::size_t i = 0; i < M; ++i) {
        if (L.stations[i].sched != SchedStrategy::FCFS) continue;
        const T mi = num_traits<T>::from_double(L.stations[i].nservers);
        const T rho = station_rho(L, i, a1, mi);
        if (!(rho < T(one - tol))) {
            for (std::size_t r = 0; r < K; ++r) {
                Q(i, r) = num_traits<T>::from_double(L.classes[r].population);
                R(i, r) = (Tp(i, r) > zero) ? T(Q(i, r) / Tp(i, r)) : zero;
            }
            continue;
        }
        const Mmap<T> arv = arrival_mmap(a1, a2, i, K, false, opt.space_max);
        std::vector<PhService<T>> sl;
        for (std::size_t r = 0; r < K; ++r) sl.push_back(svc[i][r]);

        if (std::isfinite(L.cap[i])) {
            const std::size_t capK = static_cast<std::size_t>(std::llround(L.cap[i]));
            T meanQ = zero, lossProb = zero;
            const MmckDetection<T> det = mam_detect_mmck(L, i + 1, arv);  // 1-based station index
            if (det.isMmck) {
                T lamTot = zero;
                for (std::size_t r = 0; r < K; ++r)
                    if (!L.disabled[i][r]) lamTot += a1(i, r);
                const qsys::MmckResult<T> ex = qsys::qsys_mmck(
                    lamTot, det.muRate, static_cast<unsigned>(std::llround(L.stations[i].nservers)),
                    static_cast<unsigned>(capK));
                meanQ = ex.meanQueueLength;
                lossProb = ex.lossProbability;
            } else {
                const basic_detail::TruncRenorm<T> tr =
                    basic_detail::truncate_renorm(arv, sl, capK);
                meanQ = tr.meanQ;
                lossProb = tr.lossProb;
            }
            // Under FCFS the wait in queue is common to every class, so the
            // aggregate mean queue length yields one Wq and R_k = Wq + S_k.
            std::vector<T> eff(K, zero);
            T sumT = zero;
            for (std::size_t r = 0; r < K; ++r) {
                const T inflow = L.disabled[i][r] ? zero : a1(i, r);
                eff[r] = T(inflow * T(one - lossProb));
                sumT += eff[r];
            }
            T Wq = zero;
            if (sumT > zero) {
                T sw = zero;
                for (std::size_t r = 0; r < K; ++r)
                    if (!L.disabled[i][r]) sw += T(eff[r] * S(i, r));
                const T w = T(T(meanQ / sumT) - T(sw / sumT));
                Wq = (w > zero) ? w : zero;
            }
            for (std::size_t r = 0; r < K; ++r) {
                Tp(i, r) = eff[r];
                U(i, r) = T(Tp(i, r) * S(i, r) / mi);
                if (Tp(i, r) > zero) {
                    R(i, r) = T(Wq + S(i, r));
                    Q(i, r) = T(Tp(i, r) * R(i, r));
                } else {
                    R(i, r) = zero;
                    Q(i, r) = zero;
                }
            }
        } else {
            const std::vector<T> m = mmapph1fcfs_ncmean(arv, sl);
            for (std::size_t r = 0; r < K; ++r) {
                Q(i, r) = m[arv.classes() == 1 ? 0 : r];
                R(i, r) = (Tp(i, r) > zero) ? T(Q(i, r) / Tp(i, r)) : zero;
            }
        }
    }

    mva::MvaSolution<T> out;
    out.Q = Q;
    out.U = U;
    out.R = R;
    out.Tp = Tp;
    out.C.assign(K, zero);
    // X is left at zero: the reference never assigns it. See the file header.
    out.X.assign(K, zero);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < M; ++i) out.C[k] += R(i, k);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (out.Q(i, k) < zero) out.Q(i, k) = T(-out.Q(i, k));
    zero_nans(out.Q);
    zero_nans(out.U);
    zero_nans(out.R);
    for (std::size_t k = 0; k < K; ++k)
        if (std::isnan(num_traits<T>::to_double(out.C[k]))) out.C[k] = zero;
    out.method = "mna";
    out.iter = static_cast<int>(fr.iterations);
    out.lG = 0.0;
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of `solver_mna_closed.m`.
 *
 * @param L   the refreshed struct; closed chains only
 * @param opt SolverMAM's options; `tol` drives the saturation test and doubles
 *            as both fixed points' `iter_tol`
 */
template <class T>
mva::MvaSolution<T> solver_mna_closed(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)opt;
        throw UnsupportedError(
            "solver_mna_closed: the throughput bisection and the inner flow fixed point stop on a "
            "tolerance, the APH arrival fit takes a ceiling of a real reciprocal, and the "
            "MMAP[K]/PH[K]/1 station solve runs the ADDA doubling iteration; rerun this model "
            "with --arith double or --arith real");
    } else {
    using namespace mna_detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;
    const T tol = num_traits<T>::from_double(opt.tol);
    // The reference's local `config.space_max = 16`, the order budget the
    // per-station arrival superposition is held under.
    const std::size_t space_max = 16;

    check_gates(L, "solver_mna_closed");
    if (C != K)
        throw UnsupportedError(
            "solver_mna_closed: the reference drives its bisection over classes but stores the "
            "throughput in the chain-indexed lambda, and renormalizes chain c's queue lengths "
            "with the class-indexed sn.njobs(c); both are only correct when each chain holds "
            "exactly one class, and this model has " +
            std::to_string(C) + " chains over " + std::to_string(K) + " classes");
    std::vector<double> Npop(K, 0.0);
    for (std::size_t r = 0; r < K; ++r) {
        if (!std::isfinite(L.classes[r].population))
            throw UnsupportedError(
                "solver_mna_closed: MNA's closed analyzer brackets each class's throughput by its "
                "population; an open class has none and belongs to solver_mna_open");
        Npop[r] = L.classes[r].population;
    }

    Matrix<T> S(M, K, zero), scv(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.disabled[i][r] && L.rates(i, r) > zero) S(i, r) = T(one / L.rates(i, r));
            const double v = num_traits<T>::to_double(L.scv(i, r));
            scv(i, r) = std::isnan(v) ? zero : L.scv(i, r);
        }
    const Matrix<T> V = station_visits(L);
    const std::vector<std::vector<PhService<T>>> svc = service_laws(L);
    std::vector<bool> is_source(M, false);
    for (std::size_t i = 0; i < M; ++i)
        is_source[i] = (L.stations[i].nodetype == qn::NodeType::Source);
    // solver_mna_closed.m applies no round-robin correction: the deterministic
    // split is carried by the open traffic equations only, and a closed model
    // that dispatches round-robin is refused before it reaches here
    // (check_model_method). The all-ones degree keeps the shared helpers on the
    // Markovian branch, which is what the reference computes.
    const Matrix<T> kRR_one(M, K, one);

    // ---- the bisection bracket -------------------------------------------
    // The upper bound is the slowest service the class can meet at a station
    // that can queue: at that rate the class saturates whatever the routing.
    std::vector<T> lambda_lb(K, zero), lambda_ub(K, zero), lambda(K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        bool any = false;
        double best = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            if (!std::isfinite(L.stations[i].nservers)) continue;
            if (L.disabled[i][r]) continue;  // MATLAB's NaN rate, dropped by min
            const double v = num_traits<T>::to_double(L.rates(i, r));
            if (!any || v < best) {
                best = v;
                any = true;
            }
        }
        if (!any)
            throw UnsupportedError(
                "solver_mna_closed: class '" + L.classes[r].name +
                "' is served at no finite-server station, so the reference's throughput upper "
                "bound min(sn.rates(sn.nservers<Inf,k)) is empty and the assignment fails");
        lambda_ub[r] = num_traits<T>::from_double(best);
    }

    Matrix<T> Q(M, K, zero), U(M, K, zero), R(M, K, zero), Tp(M, K, zero);
    Matrix<T> a1(M, K, zero), a2(M, K, zero), f2(M * K, M * K, zero);
    std::vector<T> d2(M, zero);
    std::vector<T> QN(K, zero);
    std::vector<T> QNc(K, zero);
    for (std::size_t r = 0; r < K; ++r) QNc[r] = num_traits<T>::from_double(Npop[r]);

    // The maximum queue length any single class can reach, and hence the level
    // the queue-length distribution is evaluated up to.
    std::size_t maxLevel = 1;
    for (std::size_t r = 0; r < K; ++r) maxLevel += static_cast<std::size_t>(std::llround(Npop[r]));

    // ---- the inner flow fixed point --------------------------------------
    auto flow_sweep = [&](const std::vector<T>&,
                          std::size_t itnum) -> std::pair<std::vector<T>, std::vector<T>> {
        std::vector<T> xref(2 * M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k) {
                xref[i * K + k] = a1(i, k);
                xref[M * K + i * K + k] = a2(i, k);
            }

        // THE RENORMALIZATION IS NaN ON THE FIRST SWEEP, and stays NaN at every
        // FCFS station for the whole inner loop: Q starts at zero, so this is
        // N * 0 / 0. Nothing downstream depends on it -- the outer sweep
        // overwrites every FCFS row and the INF/PS branches overwrite their own
        // -- but the arithmetic is reproduced rather than guarded, because a
        // guard would feed the FCFS branch a zero it never sees in MATLAB.
        for (std::size_t c = 0; c < C; ++c) {
            T colsum = zero;
            for (std::size_t i = 0; i < M; ++i) colsum += Q(i, c);
            for (std::size_t i = 0; i < M; ++i)
                Q(i, c) = num_traits<T>::from_double(Npop[c] * num_traits<T>::to_double(Q(i, c)) /
                                                     num_traits<T>::to_double(colsum));
        }

        if (itnum == 1)
            for (std::size_t c = 0; c < C; ++c)
                for (std::size_t m = 0; m < M; ++m)
                    for (std::size_t k : L.inchain[c]) Tp(m, k - 1) = T(V(m, k - 1) * lambda[c]);

        superpose(L, Tp, f2, a1, a2);

        for (std::size_t i = 0; i < M; ++i) {
            if (L.stations[i].nodetype == qn::NodeType::Join) continue;  // no-op in the reference
            const SchedStrategy sched = L.stations[i].sched;
            const T mi = num_traits<T>::from_double(L.stations[i].nservers);
            if (sched == SchedStrategy::INF) {
                d2[i] = a2(i, 0);
                for (std::size_t c = 0; c < C; ++c)
                    for (std::size_t k : L.inchain[c]) {
                        const std::size_t r = k - 1;
                        Tp(i, r) = a1(i, r);
                        U(i, r) = T(S(i, r) * Tp(i, r));
                        Q(i, r) = T(Tp(i, r) * S(i, r) * V(i, r));
                        R(i, r) = (Tp(i, r) > zero) ? T(Q(i, r) / Tp(i, r)) : zero;
                    }
            } else if (sched == SchedStrategy::PS) {
                using std::pow;
                for (std::size_t c = 0; c < C; ++c) {
                    double Nc = 0.0;
                    for (std::size_t k : L.inchain[c]) Nc += Npop[k - 1];
                    for (std::size_t k : L.inchain[c]) {
                        const std::size_t r = k - 1;
                        Tp(i, r) = T(lambda[c] * V(i, r));
                        U(i, r) = T(S(i, r) * Tp(i, r));
                    }
                    T usum = zero;
                    for (std::size_t r = 0; r < K; ++r) usum += U(i, r);
                    const T ftol = num_traits<T>::from_double(GlobalConstants::FineTol);
                    const T uden = (usum < T(one - ftol)) ? usum : T(one - ftol);
                    for (std::size_t k : L.inchain[c]) {
                        const std::size_t r = k - 1;
                        // The finite-population geometric bound: the U^(N+1)
                        // term is what keeps a closed class's queue length from
                        // running past its own population.
                        const T tail = pow(U(i, r), num_traits<T>::from_double(Nc + 1.0));
                        Q(i, r) = T(T(U(i, r) - tail) / T(one - uden));
                        R(i, r) = (Tp(i, r) > zero) ? T(Q(i, r) / Tp(i, r)) : zero;
                    }
                }
            } else if (sched == SchedStrategy::FCFS) {
                T lambda_ist = zero;
                for (std::size_t r = 0; r < K; ++r) lambda_ist += a1(i, r);
                const T rho = station_rho(L, i, a1, mi);
                if (rho < T(one - tol)) {
                    d2[i] = qna_departure_scv(L, i, a1, a2, scv, lambda_ist, rho, mi);
                } else {
                    for (std::size_t r = 0; r < K; ++r) Q(i, r) = num_traits<T>::from_double(Npop[r]);
                    d2[i] = one;
                }
                for (std::size_t r = 0; r < K; ++r) {
                    Tp(i, r) = a1(i, r);
                    U(i, r) = T(Tp(i, r) * S(i, r) / mi);
                    R(i, r) = (Tp(i, r) > zero) ? T(Q(i, r) / Tp(i, r)) : zero;
                }
            } else if (sched != SchedStrategy::EXT) {
                throw UnsupportedError(
                    std::string("solver_mna_closed: no isolated-station solution for ") +
                    lang::sched_to_text(sched) + " scheduling at station '" +
                    L.stations[i].name + "'");
            }
        }

        split(L, is_source, d2, kRR_one, f2);

        std::vector<T> xnew(2 * M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k) {
                xnew[i * K + k] = a1(i, k);
                xnew[M * K + i * K + k] = a2(i, k);
            }
        return std::make_pair(xnew, xref);
    };

    // ---- the outer bisection on the per-class throughput -------------------
    auto outer_sweep = [&](const std::vector<T>&,
                           std::size_t itout) -> std::pair<std::vector<T>, std::vector<T>> {
        if (itout != 1) {
            // Too few jobs in the network means the rate was too low, so the
            // current iterate becomes the new lower bracket, and conversely.
            for (std::size_t r = 0; r < K; ++r) {
                if (QN[r] < QNc[r]) lambda_lb[r] = lambda[r];
                else lambda_ub[r] = lambda[r];
                lambda[r] = T(T(lambda_ub[r] + lambda_lb[r]) / num_traits<T>::from_int(2));
            }
        } else {
            lambda = lambda_ub;
        }

        Q = Matrix<T>(M, K, zero);
        U = Matrix<T>(M, K, zero);
        R = Matrix<T>(M, K, zero);
        Tp = Matrix<T>(M, K, zero);
        a1 = Matrix<T>(M, K, zero);
        a2 = Matrix<T>(M, K, zero);
        d2.assign(M, zero);
        f2 = init_flow_scv(L, is_source, kRR_one);

        da::FpiOptions io;
        io.iter_max = static_cast<std::size_t>(opt.iter_max) + 1;
        io.iter_tol = opt.tol;
        io.nanstop = true;
        da::da_fpi<T>(flow_sweep, std::vector<T>(2 * M * K, zero), io);

        for (std::size_t i = 0; i < M; ++i) {
            if (L.stations[i].sched != SchedStrategy::FCFS) continue;
            const T mi = num_traits<T>::from_double(L.stations[i].nservers);
            const T rho = station_rho(L, i, a1, mi);
            if (!(rho < T(one - tol))) {
                for (std::size_t r = 0; r < K; ++r) Q(i, r) = num_traits<T>::from_double(Npop[r]);
            } else {
                const Mmap<T> arv = arrival_mmap(a1, a2, i, K, true, space_max);
                Map<T> probe;
                probe.D0 = arv.D0;
                probe.D1 = arv.Dc[0];
                if (num_traits<T>::to_double(map_lambda(probe)) < GlobalConstants::FineTol) {
                    for (std::size_t r = 0; r < K; ++r)
                        Q(i, r) = (L.rates(i, 0) > zero)
                                      ? T(num_traits<T>::from_double(GlobalConstants::FineTol) /
                                          L.rates(i, 0))
                                      : zero;
                } else {
                    std::vector<PhService<T>> sl;
                    for (std::size_t r = 0; r < K; ++r) sl.push_back(svc[i][r]);
                    const std::vector<std::vector<T>> pd = mmapph1fcfs_ncdistr(arv, sl, maxLevel);
                    // CLASS 1's marginal is the only one read: the reference
                    // captures a single output from an analyzer that returns
                    // one distribution per class, and truncates that same
                    // vector at each class's own population.
                    const std::vector<T>& pdistr = pd[0];
                    T head = zero;
                    for (std::size_t n = 0; n + 1 < maxLevel; ++n) head += pdistr[n];
                    for (std::size_t r = 0; r < K; ++r) {
                        const std::size_t Nk = static_cast<std::size_t>(std::llround(Npop[r]));
                        std::vector<T> p(Nk + 1, zero);
                        for (std::size_t n = 0; n <= Nk; ++n) p[n] = num_abs(pdistr[n]);
                        p[Nk] = num_abs(T(one - head));
                        T mass = zero;
                        for (std::size_t n = 0; n <= Nk; ++n) mass += p[n];
                        T m = zero;
                        if (mass > zero)
                            for (std::size_t n = 0; n <= Nk; ++n)
                                m += num_traits<T>::from_int((int)n) * T(p[n] / mass);
                        if (m < zero) m = zero;
                        if (num_traits<T>::to_double(m) > static_cast<double>(Nk))
                            m = num_traits<T>::from_int((int)Nk);
                        Q(i, r) = m;
                    }
                }
            }
            for (std::size_t r = 0; r < K; ++r)
                R(i, r) = (Tp(i, r) > zero) ? T(Q(i, r) / Tp(i, r)) : zero;
        }

        for (std::size_t r = 0; r < K; ++r) {
            QN[r] = zero;
            for (std::size_t i = 0; i < M; ++i) QN[r] += Q(i, r);
        }
        return std::make_pair(QN, QNc);
    };

    da::FpiOptions fo;
    fo.iter_max = static_cast<std::size_t>(opt.iter_max);
    fo.iter_tol = opt.tol;
    fo.nanstop = true;
    const da::FpiResult<T> fr = da::da_fpi<T>(outer_sweep, std::vector<T>(K, zero), fo);

    // ---- the terminal renormalization ------------------------------------
    for (std::size_t c = 0; c < C; ++c) {
        T colsum = zero;
        for (std::size_t i = 0; i < M; ++i) colsum += Q(i, c);
        for (std::size_t i = 0; i < M; ++i)
            Q(i, c) = num_traits<T>::from_double(Npop[c] * num_traits<T>::to_double(Q(i, c)) /
                                                 num_traits<T>::to_double(colsum));
    }
    // An infinite server's utilization IS its queue length.
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].sched == SchedStrategy::INF)
            for (std::size_t r = 0; r < K; ++r) U(i, r) = Q(i, r);

    mva::MvaSolution<T> out;
    out.Q = Q;
    out.U = U;
    out.R = R;
    out.Tp = Tp;
    out.C.assign(K, zero);
    // X is left at zero: the reference never assigns it. See the file header.
    out.X.assign(K, zero);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < M; ++i) out.C[k] += R(i, k);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (out.Q(i, k) < zero) out.Q(i, k) = T(-out.Q(i, k));
    zero_nans(out.Q);
    zero_nans(out.U);
    zero_nans(out.R);
    for (std::size_t k = 0; k < K; ++k)
        if (std::isnan(num_traits<T>::to_double(out.C[k]))) out.C[k] = zero;
    out.method = "mna";
    out.iter = static_cast<int>(fr.iterations);
    out.lG = 0.0;
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MNA_H
