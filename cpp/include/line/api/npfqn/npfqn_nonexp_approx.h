/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_NONEXP_APPROX_H
#define LINE_API_NPFQN_NONEXP_APPROX_H

/**
 * Handler for non-exponential service and arrival processes in AMVA and NC.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_nonexp_approx.m, cross-checked
 * against jar/src/main/java/jline/api/npfqn/Npfqn_nonexp_approx.java. See the
 * note on the method names at the end of this comment: the two disagree there.
 *
 * The 'interp' method replaces the service time of every FCFS station whose
 * per-class demands or SCVs make it non-product-form by the WSC 2020
 * interpolation (LINE paper, Sec. 4.2) between the M/G/1 diffusion decay rate
 *
 *   eta_i = exp(-2 (1 - rho_i) / (c2_{s,i} + c2_{a,i} rho_i))     (Kobayashi)
 *
 * and the multiserver asymptotic decay rate gamma_i = (rho_i^{c_i} + rho_i)/2,
 * with weights a_i = b_i = rho_i^8. The multiserver effect is absorbed into
 * the scaled service time, so the station is returned with one server.
 *
 * Arithmetic. eta carries an exp, so this requires transcendental arithmetic
 * and cannot be instantiated at T = Rational. rho^{c_i} is written with a real
 * exponent because MATLAB stores the server count as a double; rho^8 uses the
 * integer power, which is the same value in every arithmetic.
 *
 * Model layer. The MATLAB function reads exactly two fields of the
 * NetworkStruct, sn.sched(i) (only ever compared against SchedStrategy.FCFS)
 * and sn.rates(i,k) (only ever compared against zero). Since the model layer
 * is not part of this port, those two are passed explicitly as isFCFS and
 * rates; nothing else about sn is consulted by the algorithm, so this is a
 * transcription of the same code, not a reduced variant.
 *
 * MATLAB / JAR disagreement (unresolved here, reported upstream): MATLAB's
 * no-op branch is case {'default','none'} plus case {'hvmva'}, and
 * solver_amvald_forward.m dispatches on the string 'hvmva'. The JAR instead
 * accepts "hmva" and throws IllegalArgumentException on anything else, while
 * its own Solver_amvald.java tests options.config.highvar against "hvmva", so
 * a JAR run with highvar='hvmva' reaching SolverFluid throws instead of
 * no-opping. Python native (line_solver/api/npfqn/nonexp.py) copies the JAR.
 * This port follows MATLAB, the reference implementation, and accepts both
 * spellings so that a caller ported from either side behaves identically.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/npfqn/npfqn_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

/** Return value, mirroring MATLAB's [ST,gamma,nservers,rho,scva,scvs,eta]. */
template <class T>
struct NonexpApproxResult {
    Matrix<T> ST;             ///< (M x R) scaled service times
    std::vector<T> gamma;     ///< (M) multiserver asymptotic decay rate
    std::vector<T> nservers;  ///< (M) server counts, set to 1 where rescaled
    std::vector<T> rho;       ///< (M) station utilization
    std::vector<T> scva;      ///< (M) arrival SCV used by the approximation
    std::vector<T> scvs;      ///< (M) throughput-weighted service SCV
    std::vector<T> eta;       ///< (M) diffusion decay rate
};

/**
 * @param method   'default', 'none', 'hvmva' (all no-ops) or 'interp'
 * @param isFCFS   (M) sn.sched(i) == SchedStrategy.FCFS
 * @param rates    (M x R) sn.rates
 * @param ST       (M x R) service times
 * @param V        (M x R) visit ratios; present in the MATLAB signature but
 *                 never read by it, kept here for a 1:1 argument list
 * @param SCV      (M x R) service SCVs
 * @param Tput     (M x R) per-class throughputs
 * @param U        (M x R) per-class utilizations
 * @param gamma    (M) input decay rates
 * @param nservers (M) input server counts
 */
template <class T>
NonexpApproxResult<T> npfqn_nonexp_approx(const std::string& method, const std::vector<bool>& isFCFS,
                                          const Matrix<T>& rates, const Matrix<T>& ST,
                                          const Matrix<T>& V, const Matrix<T>& SCV,
                                          const Matrix<T>& Tput, const Matrix<T>& U,
                                          const std::vector<T>& gamma,
                                          const std::vector<T>& nservers) {
    static_assert(num_traits<T>::has_transcendental,
                  "npfqn_nonexp_approx requires transcendental arithmetic");
    (void)V;
    const std::size_t M = isFCFS.size();
    const std::size_t R = ST.cols();
    if (ST.rows() != M || SCV.rows() != M || Tput.rows() != M || U.rows() != M || rates.rows() != M)
        throw InputError("npfqn_nonexp_approx: an input matrix has the wrong number of stations");
    if (gamma.size() != M || nservers.size() != M)
        throw InputError("npfqn_nonexp_approx: gamma or nservers has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T fineTol = num_traits<T>::from_double(1e-8);   // GlobalConstants.FineTol
    const T unitTol = num_traits<T>::from_double(1e-6);

    NonexpApproxResult<T> res;
    res.ST = ST;
    res.gamma = gamma;
    res.nservers = nservers;
    res.rho.assign(M, zero);
    res.scva.assign(M, one);
    res.scvs.assign(M, one);
    res.eta.assign(M, one);

    if (method == "default" || method == "none" || method == "hvmva" || method == "hmva") return res;
    if (method != "interp")
        throw InputError("npfqn_nonexp_approx: unsupported method '" + method + "'");

    for (std::size_t ist = 0; ist < M; ++ist) {
        std::vector<std::size_t> nnz;
        for (std::size_t k = 0; k < R; ++k)
            if (detail::num_isfinite(res.ST(ist, k)) && detail::num_isfinite(SCV(ist, k)))
                nnz.push_back(k);
        for (std::size_t k : nnz) res.rho[ist] += U(ist, k);
        if (nnz.empty() || !isFCFS[ist]) continue;

        // non-product-form test: unequal demands, or any SCV away from 1
        T stMin = res.ST(ist, nnz[0]), stMax = res.ST(ist, nnz[0]);
        T scvMin = SCV(ist, nnz[0]), scvMax = SCV(ist, nnz[0]);
        for (std::size_t k : nnz) {
            if (res.ST(ist, k) < stMin) stMin = res.ST(ist, k);
            if (res.ST(ist, k) > stMax) stMax = res.ST(ist, k);
            if (SCV(ist, k) < scvMin) scvMin = SCV(ist, k);
            if (SCV(ist, k) > scvMax) scvMax = SCV(ist, k);
        }
        const bool nonPf = stMax - stMin > zero || scvMax > one + fineTol || scvMin < one - fineTol;
        if (!nonPf) continue;

        res.scva[ist] = one;  // use an M/G/k approximation
        T tsum = zero, wsum = zero;
        for (std::size_t k : nnz) {
            tsum += Tput(ist, k);
            wsum += SCV(ist, k) * Tput(ist, k);
        }
        res.scvs[ist] = wsum / tsum;
        // multi-server asymptotic decay rate
        res.gamma[ist] = (detail::num_pow(res.rho[ist], res.nservers[ist]) + res.rho[ist]) / two;

        if (res.scvs[ist] > one - unitTol && res.scvs[ist] < one + unitTol && res.nservers[ist] == one) {
            res.eta[ist] = res.rho[ist];  // M/M/1
        } else {
            // single-server diffusion approximation (Kobayashi, JACM)
            res.eta[ist] = detail::num_exp(
                T(-two * (one - res.rho[ist]) / (res.scvs[ist] + res.scva[ist] * res.rho[ist])));
        }

        // interpolation (Sec. 4.2, LINE paper at WSC 2020). The ai, bi
        // coefficients use the 8th power, which numerically beats the 4th.
        const T ai = num_pow_int(res.rho[ist], 8u);
        const T bi = ai;
        T oneMinusAi = one - ai;
        if (oneMinusAi < zero) oneMinusAi = zero;
        T oneMinusBi = one - bi;
        if (oneMinusBi < zero) oneMinusBi = zero;
        for (std::size_t k : nnz) {
            if (rates(ist, k) > zero)
                res.ST(ist, k) = oneMinusAi * res.ST(ist, k) +
                                 ai * (bi * res.eta[ist] + oneMinusBi * res.gamma[ist]) *
                                     (res.nservers[ist] / tsum);
        }
        // multi-server effects are now inside the scaled service times
        res.nservers[ist] = one;
    }
    return res;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_NONEXP_APPROX_H
