/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LINEARIZERMX_H
#define LINE_API_PFQN_LINEARIZERMX_H

/**
 * Linearizer for mixed open/closed queueing networks.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_linearizermx.m, cross-checked
 * against jar/src/main/java/jline/api/pfqn/mva/Pfqn_linearizermx.java.
 *
 * The open classes are solved in closed form from their fixed arrival rates:
 * X_r = lambda_r and U(i,r) = lambda_r L(i,r). Their aggregate utilization
 * U^o(i) = sum_{r open} U(i,r) then inflates the closed-class demands by the
 * Bard-Schweitzer/Reiser open-class correction
 *
 *   Dc(i,c) = L(i,c) / (1 - U^o(i)),
 *
 * and the resulting purely closed subnetwork is handed to one of the
 * Linearizer variants: pfqn_linearizer, pfqn_gflinearizer with alpha = 2, or
 * pfqn_egflinearizer with a Gompertz alpha per closed class, or
 * pfqn_linearizerms when any station has more than one server. The open-class
 * residence times are finally recovered from the closed-class queue lengths as
 * W(i,r) = L(i,r) (1 + sum_c Qc(i,c)) / (1 - U^o(i)).
 *
 * Arithmetic: TRANSCENDENTAL-GATED, on two independent grounds. Every branch
 * delegates to a Linearizer variant that stops on a tolerance, so the answer
 * is a fixed point only to within tol; and the 'egflin' branch evaluates the
 * Gompertz exponent 0.6 + 1.4 exp(-8 exp(-0.8 N_c)), a genuine transcendental
 * of the closed population, which has no meaning in an exact field.
 *
 * Open classes. MATLAB marks them with N_r = Inf. There is no infinity in a
 * general ordered field, so this port marks them with a negative entry
 * (kOpenClass) in the otherwise non-negative population vector, and the
 * `lambda(isnan) = 0` style scrubbing of the reference is dropped: NaN is a
 * floating-point artefact, not a value of T, and silently rewriting a caller's
 * input to zero would hide a modelling error rather than fix one.
 *
 * MATLAB-vs-JAR disagreement on the 'egflin' alpha, resolved in the JAR's
 * favour. MATLAB builds the exponent vector in the GLOBAL class index space
 *
 *     alphaM = zeros(1,R);                       % R = total class count
 *     for ridx = 1:length(closedClasses)
 *         r = closedClasses(ridx);
 *         alphaM(r) = 0.6 + 1.4*exp(-8*exp(-0.8*N(r)));
 *     end
 *
 * but pfqn_egflinearizer consumes it in the CLOSED-class index space, since it
 * is called with Dc, which has only length(closedClasses) columns and whose
 * class r is closedClasses(r). Whenever an open class precedes a closed one
 * the two index spaces differ and the closed class silently receives the
 * leading zero of alphaM, i.e. alpha = 0, so N_c^alpha_c collapses to 1 and
 * the Gompertz scaling is switched off entirely. Verified on the two-station,
 * two-class model of cpp/tests/test_pfqn_linearizer.cpp with class 1 open at
 * lambda = 0.4 and class 2 closed at N = 3: MATLAB pfqn_linearizermx returns
 * X_2 = 0.85177305755470, which reproduces exactly a direct
 * pfqn_egflinearizer call with alpha = 0, whereas the intended exponent
 * alpha = 1.2775503648739 gives X_2 = 0.84260134365230. The JAR indexes
 * alphaM over Nclosed and is correct; this port follows the JAR. The defect is
 * invisible in the single-class-per-index case where closedClasses == 1:R.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_egflinearizer.h"
#include "line/api/pfqn/pfqn_gflinearizer.h"
#include "line/api/pfqn/pfqn_linearizer.h"
#include "line/api/pfqn/pfqn_linearizerms.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Population sentinel marking an open class, standing in for MATLAB's Inf. */
constexpr int kOpenClass = -1;

/** Which Linearizer variant solves the closed subnetwork. */
enum class LinearizerMxMethod { Lin, Gflin, Egflin };

/**
 * @param lambda   (R) per-class arrival rate; must be zero on closed classes
 * @param L        (M x R) service demands
 * @param N        (R) population per class, kOpenClass for an open class
 * @param Z        (K x R) think times, summed over rows; may be empty
 * @param nservers (M) servers per station; all-one selects the single-server
 *                 branch, anything larger routes to pfqn_linearizerms
 * @param type     (M) scheduling discipline; empty means all-PS, as in MATLAB
 * @param tol      convergence tolerance
 * @param maxiter  total inner-iteration budget
 * @param method   Linearizer variant for the closed subnetwork
 * @param QN0      warm start, (M x R) or (M x closed count); may be empty
 */
template <class T>
LinearizerResult<T> pfqn_linearizermx(const std::vector<T>& lambda, const Matrix<T>& L,
                                      const std::vector<int>& N, const Matrix<T>& Z,
                                      const std::vector<int>& nservers,
                                      const std::vector<SchedStrategy>& type, double tol,
                                      int maxiter, LinearizerMxMethod method,
                                      const Matrix<T>& QN0) {
    // runtime gating rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)

    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError(
            "pfqn_linearizermx: demand matrix and population vector disagree on the class count");
    if (lambda.size() != R)
        throw InputError("pfqn_linearizermx: arrival-rate vector has the wrong class count");
    if (nservers.size() != M)
        throw InputError("pfqn_linearizermx: server-count vector has the wrong station count");
    if (!type.empty() && type.size() != M)
        throw InputError("pfqn_linearizermx: scheduling vector has the wrong station count");
    for (int c : nservers)
        if (c < 1) throw InputError("pfqn_linearizermx: server count below one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<std::size_t> openC, closedC;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] == kOpenClass) {
            openC.push_back(r);
        } else {
            if (N[r] < 0) throw InputError("pfqn_linearizermx: negative population");
            // The reference refuses an arrival rate on a class that also has a
            // finite positive population: it is neither open nor closed.
            if (N[r] > 0 && lambda[r] != zero)
                throw InputError(
                    "pfqn_linearizermx: arrival rate cannot be specified on closed classes");
            closedC.push_back(r);
        }
    }

    LinearizerResult<T> res;
    res.Q = Matrix<T>(M, R, zero);
    res.U = Matrix<T>(M, R, zero);
    res.W = Matrix<T>(M, R, zero);
    res.C.assign(R, zero);
    res.X.assign(R, zero);
    res.totiter = 0;

    for (std::size_t r : openC) {
        res.X[r] = lambda[r];
        for (std::size_t i = 0; i < M; ++i) res.U(i, r) = lambda[r] * L(i, r);
    }
    // Aggregate open-class utilization, before any closed class contributes.
    std::vector<T> Ut(M, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r : openC) Ut[i] += res.U(i, r);

    const std::vector<T> Zs = sum_rows(Z, R);

    const std::size_t Rc = closedC.size();
    Matrix<T> Dc(M, Rc, zero);
    for (std::size_t i = 0; i < M; ++i) {
        const T slack = one - Ut[i];
        // saturated-station rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        if (!(slack > zero))
            throw NumericError(
                "pfqn_linearizermx: open-class traffic saturates a station, the closed-class "
                "demand correction 1/(1 - U_open) is not positive");
        for (std::size_t k = 0; k < Rc; ++k) Dc(i, k) = L(i, closedC[k]) / slack;
    }

    std::vector<int> Nc(Rc, 0);
    Matrix<T> Zc(1, Rc, zero);
    for (std::size_t k = 0; k < Rc; ++k) {
        Nc[k] = N[closedC[k]];
        Zc(0, k) = Zs[closedC[k]];
    }

    // Warm start, accepted either over all classes or over the closed ones.
    Matrix<T> QN0c;
    if (!QN0.empty() && QN0.rows() == M) {
        if (QN0.cols() == R) {
            QN0c = Matrix<T>(M, Rc, zero);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t k = 0; k < Rc; ++k) QN0c(i, k) = QN0(i, closedC[k]);
        } else if (QN0.cols() == Rc) {
            QN0c = QN0;
        }
    }

    int cmax = 1;
    for (int c : nservers)
        if (c > cmax) cmax = c;

    LinearizerResult<T> sub;
    if (Rc == 0) {
        sub.Q = Matrix<T>(M, 0, zero);
    } else if (cmax == 1) {
        switch (method) {
            case LinearizerMxMethod::Gflin:
                sub = pfqn_gflinearizer(Dc, Nc, Zc, type, tol, maxiter,
                                        num_traits<T>::from_double(2.0), QN0c);
                break;
            case LinearizerMxMethod::Egflin: {
                // Gompertz exponent rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                std::vector<T> alphaM(Rc, zero);
                for (std::size_t k = 0; k < Rc; ++k) {
                    if constexpr (num_traits<T>::has_transcendental) {
                        using std::exp;
                        const T n = num_traits<T>::from_int(Nc[k]);
                        const T inner = exp(T(num_traits<T>::from_double(-0.8) * n));
                        const T outer = exp(T(num_traits<T>::from_double(-8.0) * inner));
                        alphaM[k] = num_traits<T>::from_double(0.6) +
                                    num_traits<T>::from_double(1.4) * outer;
                    } else {
                        const double n = static_cast<double>(Nc[k]);
                        alphaM[k] = num_traits<T>::from_double(
                            0.6 + 1.4 * std::exp(-8.0 * std::exp(-0.8 * n)));
                    }
                }
                sub = pfqn_egflinearizer(Dc, Nc, Zc, type, tol, maxiter, alphaM, QN0c);
                break;
            }
            case LinearizerMxMethod::Lin:
            default:
                sub = pfqn_linearizer(Dc, Nc, Zc, type, tol, maxiter, QN0c);
                break;
        }
    } else {
        sub = pfqn_linearizerms(Dc, Nc, Zc, nservers, type, tol, maxiter, QN0c);
    }
    res.totiter = sub.totiter;

    for (std::size_t k = 0; k < Rc; ++k) {
        const std::size_t r = closedC[k];
        res.X[r] = sub.X[k];
        res.C[r] = sub.C[k];
        for (std::size_t i = 0; i < M; ++i) {
            res.Q(i, r) = sub.Q(i, k);
            res.W(i, r) = sub.W(i, k);
            // Recomputed from the ORIGINAL demand L, not from the corrected
            // Dc the subnetwork was solved with, exactly as in the reference.
            res.U(i, r) = res.X[r] * L(i, r);
        }
    }

    for (std::size_t i = 0; i < M; ++i) {
        T qsum = one;
        for (std::size_t k = 0; k < Rc; ++k) qsum += sub.Q(i, k);
        for (std::size_t r : openC) {
            res.W(i, r) = L(i, r) * qsum / (one - Ut[i]);
            res.Q(i, r) = res.W(i, r) * res.X[r];
        }
    }
    for (std::size_t r : openC) {
        T c = zero;
        for (std::size_t i = 0; i < M; ++i) c += res.W(i, r);
        res.C[r] = c;
    }
    return res;
}

/** MATLAB defaults: all-PS, tol = 1e-8, maxiter = 1000, 'egflin', no warm start. */
template <class T>
LinearizerResult<T> pfqn_linearizermx(const std::vector<T>& lambda, const Matrix<T>& L,
                                      const std::vector<int>& N, const Matrix<T>& Z,
                                      const std::vector<int>& nservers,
                                      LinearizerMxMethod method = LinearizerMxMethod::Egflin) {
    return pfqn_linearizermx(lambda, L, N, Z, nservers, std::vector<SchedStrategy>(), 1e-8, 1000,
                             method, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LINEARIZERMX_H
