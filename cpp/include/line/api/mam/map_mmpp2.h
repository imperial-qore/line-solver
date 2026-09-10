/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_MMPP2_H
#define LINE_API_MAM_MAP_MMPP2_H

/**
 * Fit an MMPP(2) to a mean, an SCV, a skewness and a lag-1 autocorrelation.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_mmpp2.m. An MMPP(2) is the
 * two-state Markov-modulated Poisson process
 *
 *   D0 = [-mu00-q01, q01; q10, -mu11-q10],   D1 = diag(mu00, mu11),
 *
 * and the four rates are the closed-form inverse of its first three moments and
 * its autocorrelation decay rate G2. That inverse is Maple output: two branches,
 * the second some 18 KB of algebra in a single expression.
 *
 * THE ALGEBRA IS MACHINE-TRANSCRIBED, NOT RETYPED. It was produced by
 * `cpp/tools/matlab_expr_to_cpp.py`, which parses the MATLAB expression and
 * re-emits it, then binds the repeated radicals and denominators to `cseN`
 * temporaries (Maple repeats the same 2 KB radical dozens of times, so the
 * literal form is both unreadable and O(repeats) to evaluate). The generated
 * form was checked against the MATLAB source evaluated term by term on 40
 * random feasible inputs: the relative deviation is EXACTLY zero. Do not
 * hand-edit the generated blocks; regenerate them.
 *
 * THE FEASIBILITY GATES ARE THE REFERENCE'S, and each names what it refuses. An
 * MMPP(2) is over-dispersed (SCV >= 1) and non-negatively autocorrelated, and at
 * SCV = 1 exactly the fit is degenerate: G2 divides by (1 - 1/SCV) and every
 * rate returns NaN, so that boundary is refused by name with a pointer to the
 * Poisson process it is really asking for. ACF1 = -1 and SKEW = -1 are the
 * reference's sentinels for "give me the extreme feasible value".
 *
 * ARITHMETIC: transcendental, for the radicals.
 */

#include <cmath>

#include "line/api/mam/map_dist.h"
#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param MEAN   mean inter-arrival time
 * @param SCV_in squared coefficient of variation, which must exceed one
 * @param SKEW   skewness, or -1 for the minimum-skewness fit
 * @param ACF1   lag-1 autocorrelation, or -1 for the maximum feasible value
 */
template <class T>
Map<T> map_mmpp2(const T& MEAN, const T& SCV_in, const T& SKEW, const T& ACF1) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_mmpp2 inverts the moment equations through radicals");
    using fitdetail::num_sqrt;
    using fitdetail::pw;
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T zero = num_traits<T>::from_int(0);
    const double FEASTOLD = std::pow(10.0, -static_cast<double>(map_feastol()));
    const T FEASTOL = num_traits<T>::from_double(FEASTOLD);

    const T E1 = MEAN;
    const T E2 = T((one + SCV_in) * E1 * E1);
    T E3 = T(-(two * pw(E1, 3) - num_traits<T>::from_int(3) * E1 * E2 -
               SKEW * num_sqrt(T(pw(T(E2 - E1 * E1), 3)))));

    if (SCV_in < one - FEASTOL)
        throw InputError(
            "map_mmpp2: the SCV is infeasible, the inter-arrival times of an MMPP(2) are "
            "over-dispersed (SCV >= 1)");
    if (num_abs(T(SCV_in - one)) <= FEASTOL)
        throw InputError(
            "map_mmpp2: SCV = 1 is the Poisson boundary, where the MMPP(2) fit is degenerate: "
            "the decay rate G2 = ACF1/(1 - 1/SCV)/0.5 divides by zero and every rate comes back "
            "NaN. Use map_exponential_mean for a Poisson process");

    const T RHO0MAX = T(num_traits<T>::from_rational(1, 2) * (one - one / SCV_in));
    const bool acfAuto = (ACF1 == num_traits<T>::from_int(-1));
    if (!acfAuto) {
        if (ACF1 < -FEASTOL)
            throw InputError(
                "map_mmpp2: a negative ACF1 is infeasible, an MMPP(2) cannot be negatively "
                "autocorrelated. Pass ACF1 = -1 to request the maximum feasible autocorrelation");
        if (ACF1 > RHO0MAX + FEASTOL)
            throw InputError(
                "map_mmpp2: ACF1 exceeds the maximum lag-1 autocorrelation feasible at this SCV. "
                "Pass ACF1 = -1 to request it");
    }
    const T G2 = acfAuto
                     ? T(one - num_traits<T>::from_double(10.0 * FEASTOLD))
                     : T(ACF1 / (one - one / SCV_in) / num_traits<T>::from_rational(1, 2));

    if (SKEW == num_traits<T>::from_int(-1) && SCV_in > one)
        E3 = T((num_traits<T>::from_rational(3, 2) + num_traits<T>::from_double(0.001)) * E2 * E2 /
               E1);
    const T E3MIN = T(num_traits<T>::from_rational(3, 2) * E2 * E2 / E1);
    if (E3 < E3MIN - FEASTOL)
        throw InputError(
            "map_mmpp2: the requested skewness gives a third moment below the minimum feasible at "
            "this SCV. Pass SKEW = -1 for the minimum-skewness fit");

    // The reference recomputes SCV from the moments before the algebra.
    const T SCV = T((E2 - E1 * E1) / (E1 * E1));

    T mu00v, mu11v, q01v, q10v;
    if (num_traits<T>::to_double(G2) < 1e-6) {
        // ---- BEGIN GENERATED (matlab_expr_to_cpp.py, map_mmpp2.m G2 < 1e-6) ----
            const T cse2 = ((num_traits<T>::from_int(6) * pw(E1, 3)) * SCV);
            const T cse3 = (num_traits<T>::from_int(3) * pw(E1, 3));
            const T cse0 = (((cse2 + (cse3 * pw(SCV, 2))) + cse3) - (num_traits<T>::from_int(2) * E3));
            const T cse1 = (cse2 - E3);

            const T mu00 = (((num_traits<T>::from_int(2) * cse1) / E1) / cse0);

            const T mu11 = num_traits<T>::from_int(0);

            const T q01 = (((((num_traits<T>::from_int(9) * pw(E1, 5)) * (SCV - num_traits<T>::from_int(1))) * ((pw(SCV, 2) - (num_traits<T>::from_int(2) * SCV)) + num_traits<T>::from_int(1))) / cse1) / cse0);

            const T q10 = ((((-num_traits<T>::from_int(3)) * (SCV - num_traits<T>::from_int(1))) * pw(E1, 2)) / cse1);
        // ---- END GENERATED ----
        mu00v = mu00;
        mu11v = mu11;
        q01v = q01;
        q10v = q10;
    } else {
        // ---- BEGIN GENERATED (matlab_expr_to_cpp.py, map_mmpp2.m else branch) ----
            const T cse45 = (num_traits<T>::from_int(18) * pw(E1, 6));
            const T cse35 = (cse45 * G2);
            const T cse49 = (num_traits<T>::from_int(6) * pw(E1, 3));
            const T cse39 = (cse49 * G2);
            const T cse43 = (num_traits<T>::from_int(12) * pw(E1, 3));
            const T cse44 = ((num_traits<T>::from_int(6) * G2) * SCV);
            const T cse50 = (num_traits<T>::from_int(9) * pw(E1, 6));
            const T cse24 = (((((((((((pw(E3, 2) - ((cse43 * SCV) * E3)) + (cse39 * E3)) - ((cse44 * pw(E1, 3)) * E3)) + (((num_traits<T>::from_int(18) * G2) * pw(SCV, 3)) * pw(E1, 6))) - (cse35 * pw(SCV, 2))) + (cse50 * pw(G2, 2))) + ((num_traits<T>::from_int(36) * pw(E1, 6)) * pw(SCV, 2))) + (cse35 * SCV)) - ((cse45 * SCV) * pw(G2, 2))) + ((cse50 * pw(SCV, 2)) * pw(G2, 2))) - cse35);
            const T cse34 = (cse49 * SCV);
            const T cse48 = (num_traits<T>::from_int(3) * pw(E1, 3));
            const T cse38 = (cse48 * G2);
            const T cse42 = ((-num_traits<T>::from_int(3)) * pw(E1, 3));
            const T cse23 = (((((cse42 * G2) + (cse38 * SCV)) - cse34) + E3) + num_sqrt(cse24));
            const T cse25 = ((((cse42 * pw(SCV, 2)) - cse34) - cse48) + (num_traits<T>::from_int(2) * E3));
            const T cse11 = ((cse43 * cse23) / cse25);
            const T cse0 = (cse11 * pw(G2, 2));
            const T cse1 = (cse11 * pw(G2, 3));
            const T cse2 = ((((num_traits<T>::from_int(24) * pw(E1, 3)) * cse23) / cse25) * pw(G2, 2));
            const T cse13 = (((num_traits<T>::from_int(36) * cse23) * pw(E1, 5)) / cse25);
            const T cse3 = (cse13 * pw(G2, 2));
            const T cse4 = ((((num_traits<T>::from_int(18) * cse23) * pw(E1, 5)) / cse25) * pw(G2, 3));
            const T cse18 = (((num_traits<T>::from_int(6) * cse23) * pw(E1, 2)) / cse25);
            const T cse5 = (cse18 * pw(G2, 2));
            const T cse46 = (num_traits<T>::from_int(18) * pw(E1, 3));
            const T cse6 = (((cse46 * cse23) / cse25) * G2);
            const T cse7 = ((((num_traits<T>::from_int(27) * pw(E1, 3)) * cse23) / cse25) * G2);
            const T cse17 = (((num_traits<T>::from_int(3) * cse23) * pw(E1, 2)) / cse25);
            const T cse8 = (cse17 * G2);
            const T cse22 = (num_traits<T>::from_int(9) * cse23);
            const T cse9 = (((cse22 * pw(E1, 5)) / cse25) * G2);
            const T cse15 = ((cse48 * cse23) / cse25);
            const T cse10 = (cse15 * G2);
            const T cse12 = (((num_traits<T>::from_int(12) * cse23) * pw(E1, 2)) / cse25);
            const T cse51 = (num_traits<T>::from_int(9) * pw(E1, 3));
            const T cse14 = ((cse51 * cse23) / cse25);
            const T cse16 = ((cse49 * cse23) / cse25);
            const T cse19 = (((num_traits<T>::from_int(4) * cse23) / cse25) * E3);
            const T cse20 = ((cse23 / cse25) * E3);
            const T cse21 = ((cse23 / E1) / cse25);
            const T cse27 = (cse43 * pw(G2, 2));
            const T cse26 = (cse27 * SCV);
            const T cse28 = (cse43 * pw(G2, 3));
            const T cse29 = (cse46 * pw(G2, 2));
            const T cse47 = (num_traits<T>::from_int(18) * pw(E1, 5));
            const T cse30 = (cse47 * pw(G2, 3));
            const T cse31 = ((num_traits<T>::from_int(27) * pw(E1, 5)) * pw(G2, 2));
            const T cse32 = ((num_traits<T>::from_int(6) * pw(E1, 2)) * pw(G2, 2));
            const T cse52 = (num_traits<T>::from_int(3) * pw(E1, 2));
            const T cse41 = (cse52 * G2);
            const T cse33 = (cse41 * E3);
            const T cse36 = (cse43 * G2);
            const T cse37 = (cse47 * G2);
            const T cse40 = (cse51 * G2);

            const T mu00 = (((G2 * (((((((((((((((((-num_traits<T>::from_int(4)) * E3) * G2) + (cse19 * G2)) - cse6) - (cse6 * pw(SCV, 2))) - cse27) - (cse0 * SCV)) + ((cse11 * G2) * SCV)) + (cse36 * pw(SCV, 2))) - (cse14 * SCV)) + cse15) + cse26) + (cse14 * pw(SCV, 2))) + cse36) + cse0) - (cse15 * pw(SCV, 3)))) / ((((((((((((((((((((((((((cse28 * SCV) + ((cse48 * pw(SCV, 3)) * G2)) - cse28) + (cse29 * pw(SCV, 2))) - cse38) + (cse7 * pw(SCV, 2))) - (cse40 * pw(SCV, 2))) + cse29) - cse26) + (cse40 * SCV)) - (cse1 * SCV)) - ((cse14 * pw(SCV, 3)) * G2)) - (cse2 * pw(SCV, 2))) - (cse20 * pw(SCV, 2))) + (cse19 * pw(G2, 2))) + cse1) - cse20) + ((((num_traits<T>::from_int(2) * cse23) / cse25) * E3) * SCV)) + (cse14 * G2)) + (cse2 * SCV)) - (cse7 * SCV)) + (cse16 * SCV)) - (cse11 * pw(SCV, 2))) - cse2) + (cse16 * pw(SCV, 3))) - ((num_traits<T>::from_int(4) * E3) * pw(G2, 2)))) / E1);

            const T mu11 = cse21;

            const T q01 = ((((-num_traits<T>::from_int(3)) * pw(E1, 2)) * ((((((((((((((((((((((-num_traits<T>::from_int(6)) * cse23) * pw(E1, 2)) / cse25) * SCV) + ((cse12 * G2) * SCV)) - (cse44 * pw(E1, 2))) - cse8) + (cse21 * E3)) + cse41) + (cse18 * pw(SCV, 2))) - ((((cse22 * pw(E1, 2)) / cse25) * pw(SCV, 2)) * G2)) + (cse41 * pw(SCV, 2))) - ((((E3 * cse23) / E1) / cse25) * SCV)) - (cse5 * SCV)) + (cse32 * SCV)) + (cse17 * pw(G2, 2))) - ((((G2 * cse23) / E1) / cse25) * E3)) - (cse52 * pw(G2, 2))) + ((cse17 * pw(SCV, 2)) * pw(G2, 2))) - ((cse52 * pw(SCV, 2)) * pw(G2, 2))) + (((((G2 * SCV) * cse23) / E1) / cse25) * E3))) / ((((((((((((((((((((((((((((((-num_traits<T>::from_int(45)) * cse23) * pw(E1, 5)) / cse25) * G2) * pw(SCV, 2)) + (((num_traits<T>::from_int(18) * pw(G2, 2)) * pw(E1, 5)) * SCV)) + cse30) - (cse31 * pw(SCV, 2))) + (cse32 * E3)) - cse31) - (cse30 * SCV)) - (cse37 * SCV)) + (cse37 * pw(SCV, 2))) + cse33) - (cse33 * SCV)) + (cse21 * pw(E3, 2))) + ((cse8 * SCV) * E3)) - (cse3 * SCV)) + cse3) + (cse13 * pw(SCV, 2))) + (((((num_traits<T>::from_int(45) * cse23) * pw(E1, 5)) / cse25) * G2) * SCV)) - ((cse12 * SCV) * E3)) - (cse8 * E3)) + (cse9 * pw(SCV, 3))) + (cse3 * pw(SCV, 2))) - (cse5 * E3)) + (cse4 * SCV)) - cse4) - cse9));

            const T q10 = ((((num_traits<T>::from_int(3) * ((((((((((((((((((cse42 * cse23) / cse25) * pw(SCV, 3)) - (cse10 * pw(SCV, 2))) + (cse49 * pw(SCV, 2))) + (cse38 * pw(SCV, 2))) + (cse15 * pw(SCV, 2))) + ((cse16 * G2) * SCV)) - (E3 * SCV)) - cse34) + (cse20 * SCV)) - (cse39 * SCV)) - (cse15 * SCV)) - cse20) + cse38) - cse10) + cse15) + E3)) * pw(E1, 2)) * ((-num_traits<T>::from_int(1)) + G2)) / cse24);
        // ---- END GENERATED ----
        mu00v = mu00;
        mu11v = mu11;
        q01v = q01;
        q10v = q10;
    }

    const T rates[4] = {mu00v, mu11v, q01v, q10v};
    for (int i = 0; i < 4; ++i) {
        const double v = num_traits<T>::to_double(rates[i]);
        if (std::isnan(v) || v < -FEASTOLD)
            throw InputError(
                "map_mmpp2: the requested (MEAN, SCV, SKEW, ACF1) is not MMPP(2)-feasible; the "
                "fit gives a rate that is not a MAP");
    }
    // The reference clamps a rate that is negative only within the tolerance.
    const T a = mu00v > zero ? mu00v : zero;
    const T b = mu11v > zero ? mu11v : zero;
    const T c = q01v > zero ? q01v : zero;
    const T d = q10v > zero ? q10v : zero;

    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = T(-a - c);
    m.D0(0, 1) = c;
    m.D0(1, 0) = d;
    m.D0(1, 1) = T(-b - d);
    m.D1(0, 0) = a;
    m.D1(1, 1) = b;
    return m;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_MMPP2_H
