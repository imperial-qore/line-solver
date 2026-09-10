/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_BLOCK_H
#define LINE_API_MAM_MAP_BLOCK_H

/**
 * Fit a MAP(2) to three moments and an autocorrelation decay rate.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_block.m and map_feasblock.m.
 * The four rates of the general (non-MMPP) MAP(2)
 *
 *   D0 = diag(-mu00-mu01, -mu10-mu11),   D1 = [mu00, mu01; mu10, mu11]
 *
 * are the closed-form inverse of (E1, E2, E3, G2). That inverse is Maple output,
 * some 57 KB of algebra across the four expressions.
 *
 * THE ALGEBRA IS MACHINE-TRANSCRIBED, NOT RETYPED, by
 * `cpp/tools/matlab_expr_to_cpp.py`, which parses the MATLAB expression, binds
 * the repeated subexpressions to `cseN` temporaries, and re-emits it. Do not
 * hand-edit the generated block; regenerate it.
 *
 * IT MUST BE EVALUATED IN COMPLEX ARITHMETIC, and this is the trap. The single
 * radicand the four expressions share goes NEGATIVE on perfectly feasible
 * moment sets -- (E1, SCV, E3/E3min, G2) = (1, 2, 2, 0.3) is one -- and the
 * imaginary parts then CANCEL in the four rates. MATLAB evaluates in complex
 * arithmetic throughout and only afterwards asks whether any entry retains an
 * imaginary part above 1e-4; that residual test, not the sign of the radicand,
 * is what rejects a fit. An earlier version of this port branched on the
 * radicand's sign and sent three of three feasible test cases to the fallback
 * while MATLAB fitted all three exactly. The generated block is therefore
 * instantiated at `fitdetail::Cplx<T>` and the reference's residual test is
 * applied to the result.
 *
 * THE FALLBACK IS THE REFERENCE'S, and it drops the third moment:
 *  - SCV >= 1: a hyperexponential-shaped MAP(2) matching E1 and E2 only, with
 *    the autocorrelation carried by the switching probability p = (1 - G2)/2 and
 *    D1 = -D0 P over P = [1-p, p; p, 1-p];
 *  - SCV < 1: the exponential of mean E1. The reference's commented-out general
 *    MAP(2) branch for SCV < 1 is NOT reinstated -- it is commented out there,
 *    so reinstating it would answer a different model than every other codebase.
 *
 * THE FALLBACK IS ITSELF BOUNDED ABOVE BY SCV 3. Its first branch rate is
 * E1 (1 - sqrt((SCV-1)/2)), which vanishes at SCV = 3 and is negative above it.
 * MATLAB returns the resulting non-generator regardless (measured in R2025a: an
 * infinite diagonal at SCV = 3 with a NaN mean, a POSITIVE diagonal at SCV = 5,
 * `map_isfeasible` 0 for both). This port refuses by name instead, since handing
 * back a matrix that is not a MAP is worse than saying so.
 *
 * `map_feasblock` is the same fit behind a moment repair: an SCV at or below one
 * is raised to 1 + tol and a third moment below (3/2) E2^2 / E1 is raised to it.
 *
 * ARITHMETIC: transcendental, for the radical.
 */

#include <cmath>
#include <cstddef>

#include "line/api/mam/map_dist.h"
#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace blockdetail {

/** The four Maple rates, evaluated in complex arithmetic as MATLAB does. */
template <class T>
struct BlockRates {
    fitdetail::Cplx<T> mu00, mu10, mu01, mu11;
};

template <class T>
BlockRates<T> block_rates(const fitdetail::Cplx<T>& E1, const fitdetail::Cplx<T>& E2,
                          const fitdetail::Cplx<T>& E3, const fitdetail::Cplx<T>& G2) {
    typedef fitdetail::Cplx<T> C;
    // Local spellings the generated block uses; see the header note.
    auto CI = [](long v) { return C(num_traits<T>::from_int(v)); };
    auto CD = [](double v) { return C(num_traits<T>::from_double(v)); };
    auto cpw = [](const C& b, unsigned e) {
        C r(num_traits<T>::from_int(1));
        for (unsigned i = 0; i < e; ++i) r = r * b;
        return r;
    };

    // ---- BEGIN GENERATED (matlab_expr_to_cpp.py, map_block.m) ----
        const C cse30 = (CI(216) * cpw(E1, 7));
        const C cse20 = (cse30 * cpw(E2, 2));
        const C cse21 = ((CI(24) * cpw(E1, 5)) * cpw(E3, 3));
        const C cse24 = ((CI(6) * E1) * cpw(E3, 3));
        const C cse28 = (CI(162) * cpw(E2, 7));
        const C cse29 = (CI(567) * cpw(E1, 2));
        const C cse31 = ((CI(270) * E1) * E3);
        const C cse34 = (CI(24) * cpw(E1, 4));
        const C cse35 = (CI(18) * cpw(E1, 2));
        const C cse47 = (CI(9) * cpw(E2, 4));
        const C cse18 = ((((((cse34 * E3) - ((CI(27) * cpw(E1, 3)) * cpw(E2, 2))) - ((cse35 * E2) * E3)) + ((CI(18) * E1) * cpw(E2, 3))) + (cpw(E3, 2) * E1)) + fitdetail::cplx_sqrt(((((((((((((((((((-CI(243)) * cpw(E1, 6)) * cpw(E2, 4)) - (cse24 * cpw(E2, 2))) + cse28) + (cse47 * cpw(E3, 2))) + cse21) + ((CI(648) * cpw(E1, 4)) * cpw(E2, 5))) - (cse29 * cpw(E2, 6))) + (cse20 * E3)) + (((CI(144) * cpw(E1, 6)) * cpw(E3, 2)) * E2)) - (((CI(756) * cpw(E1, 5)) * cpw(E2, 3)) * E3)) - (((CI(270) * cpw(E1, 4)) * cpw(E3, 2)) * cpw(E2, 2))) - (((CI(12) * cpw(E1, 3)) * cpw(E3, 3)) * E2)) + (((CI(810) * cpw(E1, 3)) * E3) * cpw(E2, 4))) + (((CI(108) * cpw(E1, 2)) * cpw(E3, 2)) * cpw(E2, 3))) - (cse31 * cpw(E2, 5))) + (cpw(E3, 4) * cpw(E1, 2)))));
        const C cse19 = ((((((((((-CI(3)) * cpw(E3, 2)) * cpw(E2, 2)) + ((CI(48) * cpw(E3, 2)) * cpw(E1, 4))) + ((CI(81) * cpw(E2, 4)) * cpw(E1, 2))) - (((CI(36) * cpw(E3, 2)) * cpw(E1, 2)) * E2)) + (((CI(90) * cpw(E2, 3)) * E3) * E1)) - (((CI(126) * cpw(E2, 2)) * cpw(E1, 3)) * E3)) - (CI(54) * cpw(E2, 5))) + ((CI(2) * E1) * cpw(E3, 3)));
        const C cse0 = (((((CI(9072) * cpw(E1, 9)) * cpw(E2, 2)) * E3) * cse18) / cse19);
        const C cse1 = (((((CI(1512) * cpw(E1, 7)) * cpw(E2, 2)) * E3) * cse18) / cse19);
        const C cse2 = (((((CI(30) * cpw(E1, 2)) * cpw(E2, 4)) * E3) * cse18) / cse19);
        const C cse3 = (((((CI(12) * E1) * cpw(E3, 3)) * cse18) / cse19) * cpw(E2, 2));
        const C cse4 = (((((CI(2) * E1) * cse18) / cse19) * cpw(E3, 2)) * cpw(E2, 3));
        const C cse5 = ((((CI(3456) * cpw(E1, 10)) * cpw(E3, 2)) * cse18) / cse19);
        const C cse6 = (((((CI(378) * E1) * E3) * cse18) / cse19) * cpw(E2, 5));
        const C cse7 = ((((CI(576) * cpw(E1, 8)) * cpw(E3, 2)) * cse18) / cse19);
        const C cse50 = (CI(4) * cpw(E1, 2));
        const C cse8 = (((cse50 * cse18) / cse19) * cpw(E3, 4));
        const C cse9 = ((((CI(9) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 4));
        const C cse10 = ((((CI(9) * E1) * cpw(E2, 6)) * cse18) / cse19);
        const C cse17 = ((CI(3) * cse18) / cse19);
        const C cse16 = (cse17 * E3);
        const C cse11 = (cse16 * cpw(E2, 5));
        const C cse12 = (((CI(5832) * cpw(E1, 8)) * cse18) / cse19);
        const C cse13 = (((CI(486) * cse18) / cse19) * cpw(E2, 8));
        const C cse14 = (((CI(162) * cse18) / cse19) * cpw(E2, 7));
        const C cse15 = (((CI(972) * cpw(E1, 6)) * cse18) / cse19);
        const C cse22 = (((CI(16) * cpw(E1, 5)) * E3) * E2);
        const C cse23 = ((CI(27) * E1) * cpw(E2, 5));
        const C cse25 = (cse47 * E3);
        const C cse26 = (CI(1728) * cpw(E1, 10));
        const C cse27 = (CI(1296) * cpw(E1, 9));
        const C cse32 = (CI(648) * cpw(E1, 3));
        const C cse33 = (CI(288) * cpw(E1, 8));
        const C cse36 = (CI(16) * cpw(E1, 7));
        const C cse37 = (CI(12) * cpw(E2, 4));
        const C cse38 = (CI(36) * cpw(E1, 3));
        const C cse39 = (CI(36) * cpw(E1, 5));
        const C cse40 = (CI(84) * cpw(E1, 4));
        const C cse41 = (CI(40) * cpw(E1, 5));
        const C cse42 = (CI(72) * cpw(E1, 6));
        const C cse43 = (CI(16) * cpw(E1, 3));
        const C cse44 = (CI(32) * cpw(E1, 7));
        const C cse45 = (CI(96) * cpw(E1, 5));
        const C cse46 = (CI(72) * cpw(E1, 2));
        const C cse48 = (CI(3) * cpw(E2, 2));
        const C cse49 = (CI(3) * cpw(E2, 5));

        const C mu00 = ((((-((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((CI(6) * cpw(E3, 4)) * cpw(E2, 2)) * cse18) / cse19) + ((CI(6) * cpw(E3, 3)) * cpw(E2, 2))) + cse13) - ((((((CI(18) * E1) * cpw(E3, 3)) * cse18) / cse19) * cpw(E2, 3)) * G2)) + ((((((CI(5994) * cpw(E1, 3)) * E3) * cse18) / cse19) * cpw(E2, 5)) * G2)) - ((((((CI(810) * E1) * E3) * cse18) / cse19) * cpw(E2, 6)) * G2)) - ((((((CI(216) * cpw(E1, 5)) * cpw(E3, 3)) * cse18) / cse19) * E2) * G2)) + ((((((CI(108) * cpw(E1, 3)) * cpw(E3, 3)) * cse18) / cse19) * cpw(E2, 2)) * G2)) + ((((((CI(162) * cpw(E1, 2)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 4)) * G2)) + ((((((CI(20088) * cpw(E1, 7)) * cpw(E2, 3)) * E3) * cse18) / cse19) * G2)) - ((((((CI(7776) * cpw(E1, 8)) * cpw(E3, 2)) * cse18) / cse19) * E2) * G2)) - (cse0 * G2)) + (((((CI(16848) * cpw(E2, 4)) * E3) * cse18) / cse19) * cpw(E1, 5))) + (cse5 * G2)) + (((((CI(9504) * cpw(E1, 8)) * cpw(E3, 2)) * cse18) / cse19) * E2)) + cse0) + (((((CI(3024) * cpw(E1, 2)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 4))) - (((((CI(288) * E1) * cpw(E3, 3)) * cse18) / cse19) * cpw(E2, 3))) + ((((cse32 * E3) * cse18) / cse19) * cpw(E2, 5))) - (((((CI(2430) * E1) * E3) * cse18) / cse19) * cpw(E2, 6))) + ((((((CI(6264) * cpw(E1, 6)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 2)) * G2)) - ((((((CI(2052) * cpw(E1, 4)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 3)) * G2)) - (((((CI(828) * cpw(E1, 3)) * cpw(E3, 3)) * cse18) / cse19) * cpw(E2, 2))) - (((((CI(3645) * cse18) / cse19) * cpw(E2, 7)) * cpw(E1, 2)) * G2)) - (((((CI(3024) * cpw(E1, 6)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 2))) - (((((CI(5940) * cpw(E1, 4)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 3))) + (((((CI(3024) * cpw(E1, 5)) * cpw(E3, 3)) * cse18) / cse19) * E2)) + (((((CI(10206) * cpw(E1, 4)) * cse18) / cse19) * cpw(E2, 6)) * G2)) + ((cse12 * G2) * cpw(E2, 4))) - (((((CI(12636) * cpw(E1, 6)) * cse18) / cse19) * G2) * cpw(E2, 5))) - ((((((CI(16524) * cpw(E2, 4)) * E3) * cse18) / cse19) * cpw(E1, 5)) * G2)) - (((((CI(24624) * cpw(E1, 7)) * cpw(E2, 3)) * E3) * cse18) / cse19)) + (cse13 * G2)) - ((((CI(168) * cpw(E1, 4)) * cse18) / cse19) * cpw(E3, 4))) + ((((CI(135) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 5))) - ((((CI(1872) * cpw(E1, 7)) * cpw(E3, 3)) * cse18) / cse19)) - ((((CI(12150) * cpw(E1, 4)) * cse18) / cse19) * cpw(E2, 6))) - (cse12 * cpw(E2, 4))) - ((((CI(4) * E1) * cse18) / cse19) * cpw(E3, 5))) + ((((CI(2187) * cse18) / cse19) * cpw(E2, 7)) * cpw(E1, 2))) + ((((CI(15552) * cpw(E1, 6)) * cse18) / cse19) * cpw(E2, 5))) + (((((CI(132) * cpw(E1, 2)) * cse18) / cse19) * cpw(E3, 4)) * E2)) + (((((CI(27) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 5)) * G2)) + (((((CI(144) * cpw(E1, 7)) * cpw(E3, 3)) * cse18) / cse19) * G2)) - cse5) + (((CI(27) * E3) * cpw(E2, 5)) * G2)) - ((CI(5184) * cpw(E1, 5)) * cpw(E2, 4))) + ((CI(4536) * cpw(E1, 7)) * cpw(E2, 3))) + ((CI(576) * cpw(E1, 7)) * cpw(E3, 2))) + ((CI(81) * E3) * cpw(E2, 5))) + (cse26 * E3)) - (cse27 * cpw(E2, 2))) - ((cse26 * G2) * E3)) + ((cse27 * G2) * cpw(E2, 2))) - (((CI(5616) * cpw(E1, 8)) * E2) * E3)) - (((CI(2592) * cpw(E1, 7)) * cpw(E2, 3)) * G2)) + (((CI(5832) * cpw(E1, 6)) * cpw(E2, 2)) * E3)) - (((CI(1080) * cpw(E1, 5)) * cpw(E3, 2)) * E2)) + (((CI(1944) * cpw(E1, 5)) * G2) * cpw(E2, 4))) - (((CI(2160) * cpw(E1, 4)) * cpw(E2, 3)) * E3)) + ((cse32 * cpw(E3, 2)) * cpw(E2, 2))) - ((cse32 * cpw(E2, 5)) * G2)) - (((CI(24) * cpw(E1, 2)) * cpw(E3, 3)) * E2)) + (((CI(54) * cpw(E1, 2)) * E3) * cpw(E2, 4))) - (((CI(126) * E1) * cpw(E3, 2)) * cpw(E2, 3))) + (((CI(81) * E1) * cpw(E2, 6)) * G2)) + (cse34 * cpw(E3, 3))) + ((CI(2430) * cpw(E1, 3)) * cpw(E2, 5))) - ((CI(405) * E1) * cpw(E2, 6))) + ((((CI(3888) * cpw(E1, 8)) * E2) * G2) * E3)) - ((((CI(3456) * cpw(E1, 6)) * E3) * G2) * cpw(E2, 2))) + ((((CI(1512) * cpw(E1, 4)) * E3) * cpw(E2, 3)) * G2)) - ((((CI(324) * cpw(E1, 2)) * E3) * G2) * cpw(E2, 4)))) / (((CI(2) * E3) * E1) - cse48)) / ((((((CI(24) * E3) * cpw(E1, 3)) - ((CI(27) * cpw(E2, 2)) * cpw(E1, 2))) - (((CI(18) * E3) * E2) * E1)) + cpw(E3, 2)) + (CI(18) * cpw(E2, 3)))) / ((((((((((-CI(12)) * cpw(E1, 4)) + ((((CI(24) * cse18) / cse19) * cpw(E1, 4)) * E3)) - ((((CI(36) * cse18) / cse19) * cpw(E1, 3)) * cpw(E2, 2))) - (((((CI(18) * cse18) / cse19) * cpw(E1, 2)) * E2) * E3)) + ((CI(12) * E2) * cpw(E1, 2))) + ((((CI(2) * cse18) / cse19) * cpw(E3, 2)) * E1)) + ((((CI(27) * cse18) / cse19) * E1) * cpw(E2, 3))) - (cse16 * cpw(E2, 2))) - cse48));

        const C mu10 = ((cse17 * ((CI(2) * cpw(E1, 2)) - E2)) * (G2 - CI(1)));

        const C mu01 = ((CI(9) * ((((((((((((((((((((((((((((((((((cse36 * E3) - ((cse37 * G2) * cpw(E1, 2))) - ((CI(12) * cpw(E2, 3)) * cpw(E1, 4))) + (cse49 * G2)) - cse49) + cse4) + (cse10 * G2)) + (cse11 * G2)) + (cse37 * cpw(E1, 2))) - ((((cse38 * cpw(E2, 5)) * cse18) / cse19) * G2)) + ((((cse39 * cpw(E2, 4)) * cse18) / cse19) * G2)) + cse2) - ((((cse40 * cse18) / cse19) * E3) * cpw(E2, 3))) + ((((cse41 * cse18) / cse19) * cpw(E3, 2)) * E2)) + (((cse38 * cse18) / cse19) * cpw(E2, 5))) - cse10) + ((((cse42 * cse18) / cse19) * E3) * cpw(E2, 2))) - ((((cse43 * cpw(E3, 2)) * cpw(E2, 2)) * cse18) / cse19)) + ((((cse44 * cpw(E3, 2)) * cse18) / cse19) * G2)) - (cse4 * G2)) - (cse2 * G2)) - (((cse44 * cse18) / cse19) * cpw(E3, 2))) + (((((cse43 * cpw(E2, 2)) * cpw(E3, 2)) * cse18) / cse19) * G2)) - (((((cse41 * E2) * cpw(E3, 2)) * cse18) / cse19) * G2)) - (((cse39 * cse18) / cse19) * cpw(E2, 4))) - cse11) + (((((cse40 * cpw(E2, 3)) * E3) * cse18) / cse19) * G2)) - (((((cse42 * cpw(E2, 2)) * E3) * cse18) / cse19) * G2)) + (((CI(4) * cpw(E2, 2)) * cpw(E1, 3)) * E3)) - ((cse36 * G2) * E3)) - cse22) + (((CI(12) * G2) * cpw(E2, 3)) * cpw(E1, 4))) + (cse22 * G2)) - ((((CI(4) * cpw(E1, 3)) * G2) * E3) * cpw(E2, 2)))) / (((((((((((((((((((((-CI(48)) * cpw(E1, 5)) * cpw(E3, 2)) + (((cse45 * cpw(E3, 3)) * cse18) / cse19)) + (((CI(108) * cpw(E1, 4)) * E3) * cpw(E2, 2))) - (((((CI(396) * cpw(E1, 4)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 2))) + ((cse38 * cpw(E3, 2)) * E2)) - (((((CI(72) * cpw(E1, 3)) * cpw(E3, 3)) * cse18) / cse19) * E2)) + (((((CI(540) * cpw(E2, 4)) * E3) * cse18) / cse19) * cpw(E1, 3))) - ((CI(54) * cpw(E1, 3)) * cpw(E2, 4))) - ((cse46 * cpw(E2, 3)) * E3)) - ((((CI(243) * cpw(E1, 2)) * cse18) / cse19) * cpw(E2, 6))) + (((((CI(288) * cpw(E1, 2)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 3))) + cse8) - (((CI(6) * E1) * cpw(E3, 2)) * cpw(E2, 2))) - cse6) - cse3) + cse23) + cse9) + cse14) + cse25));

        const C mu11 = ((-((((((((((((((((((((((((((((((((((((((((((((((((((((((((-CI(60)) * cpw(E1, 3)) * cpw(E3, 2)) * E2) + (((CI(198) * cpw(E1, 2)) * cpw(E2, 3)) * E3)) - (((CI(288) * cpw(E1, 4)) * E3) * cpw(E2, 2))) + cse25) - ((cse47 * G2) * E3)) - ((((CI(324) * cpw(E1, 4)) * G2) * E3) * cpw(E2, 2))) + ((((CI(90) * cpw(E1, 2)) * G2) * cpw(E2, 3)) * E3)) - ((CI(189) * E1) * cpw(E2, 5))) - (((CI(12) * E1) * cpw(E3, 2)) * cpw(E2, 2))) - ((cse33 * G2) * E3)) + ((cse30 * G2) * cpw(E2, 2))) - (((CI(216) * cpw(E1, 6)) * E3) * E2)) - (((CI(324) * cpw(E1, 5)) * cpw(E2, 3)) * G2)) + (((CI(162) * cpw(E1, 3)) * G2) * cpw(E2, 4))) - (cse23 * G2)) - ((((CI(120) * cpw(E1, 5)) * cpw(E3, 3)) * cse18) / cse19)) - cse14) + (((((CI(84) * cpw(E1, 3)) * cpw(E3, 3)) * cse18) / cse19) * E2)) - ((((CI(81) * cpw(E1, 2)) * cse18) / cse19) * cpw(E2, 6))) - cse8) - cse9) - cse7) - (cse15 * cpw(E2, 4))) + (((((CI(216) * cpw(E1, 4)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 2))) - (((((CI(306) * cpw(E1, 2)) * cpw(E3, 2)) * cse18) / cse19) * cpw(E2, 3))) + cse6) + cse3) + ((((CI(1134) * cpw(E1, 4)) * cse18) / cse19) * cpw(E2, 5))) - (((cse28 * cse18) / cse19) * G2)) - ((((((CI(1458) * cpw(E1, 3)) * E3) * cse18) / cse19) * G2) * cpw(E2, 4))) - (((((cse46 * cpw(E3, 2)) * cse18) / cse19) * G2) * cpw(E2, 3))) + ((((cse24 * cse18) / cse19) * G2) * cpw(E2, 2))) + ((cse15 * G2) * cpw(E2, 4))) + (((cse21 * cse18) / cse19) * G2)) - (((((CI(1836) * cpw(E1, 5)) * cpw(E2, 3)) * E3) * cse18) / cse19)) - (((((CI(1620) * cpw(E1, 4)) * cse18) / cse19) * G2) * cpw(E2, 5))) + (((((CI(891) * cpw(E1, 2)) * cse18) / cse19) * cpw(E2, 6)) * G2)) - (((((CI(9) * cse18) / cse19) * cpw(E3, 2)) * cpw(E2, 4)) * G2)) + (cse7 * G2)) + ((((((CI(540) * cpw(E1, 4)) * cpw(E3, 2)) * cse18) / cse19) * G2) * cpw(E2, 2))) - ((((((CI(24) * cpw(E1, 3)) * cpw(E3, 3)) * cse18) / cse19) * G2) * E2)) + ((((cse31 * cse18) / cse19) * G2) * cpw(E2, 5))) + cse1) + (((((CI(720) * cpw(E1, 6)) * E2) * cpw(E3, 2)) * cse18) / cse19)) + ((((((CI(2592) * cpw(E1, 5)) * cpw(E2, 3)) * E3) * cse18) / cse19) * G2)) - (cse1 * G2)) - ((((((CI(1008) * cpw(E1, 6)) * E2) * cpw(E3, 2)) * cse18) / cse19) * G2)) + ((((CI(504) * cpw(E1, 6)) * E3) * E2) * G2)) + (cse45 * cpw(E3, 2))) + ((CI(378) * cpw(E1, 3)) * cpw(E2, 4))) + (cse33 * E3)) - cse20) + (cse50 * cpw(E3, 3)))) / ((((((((((((((((((((((((-CI(96)) * cpw(E3, 2)) * cpw(E1, 6)) + ((((CI(192) * cpw(E1, 6)) * cpw(E3, 3)) * cse18) / cse19)) - (((((CI(648) * cpw(E1, 5)) * cse18) / cse19) * cpw(E3, 2)) * cpw(E2, 2))) + (((CI(108) * E3) * cpw(E2, 2)) * cpw(E1, 5))) + (((((CI(702) * cpw(E1, 4)) * cse18) / cse19) * E3) * cpw(E2, 4))) - (((((CI(192) * cpw(E1, 4)) * cpw(E3, 3)) * cse18) / cse19) * E2)) + (((CI(72) * cpw(E3, 2)) * cpw(E1, 4)) * E2)) - ((CI(4) * cpw(E3, 3)) * cpw(E1, 3))) + ((((CI(8) * cpw(E1, 3)) * cse18) / cse19) * cpw(E3, 4))) - ((((CI(243) * cpw(E1, 3)) * cse18) / cse19) * cpw(E2, 6))) + (((((CI(594) * cpw(E1, 3)) * cse18) / cse19) * cpw(E3, 2)) * cpw(E2, 3))) + ((((cse35 * cpw(E3, 3)) * cse18) / cse19) * cpw(E2, 2))) - ((CI(81) * cpw(E2, 5)) * cpw(E1, 2))) - ((((cse29 * cpw(E2, 5)) * E3) * cse18) / cse19)) - (((((CI(2) * E1) * cpw(E3, 4)) * E2) * cse18) / cse19)) - (((((CI(81) * E1) * cse18) / cse19) * cpw(E3, 2)) * cpw(E2, 4))) + ((((CI(162) * E1) * cpw(E2, 7)) * cse18) / cse19)) - (((CI(54) * E3) * cpw(E2, 4)) * E1)) + ((((CI(54) * cse18) / cse19) * E3) * cpw(E2, 6))) + (CI(54) * cpw(E2, 6))) + ((CI(3) * cpw(E2, 3)) * cpw(E3, 2))) + ((((CI(3) * cpw(E3, 3)) * cse18) / cse19) * cpw(E2, 3))));
    // ---- END GENERATED ----

    BlockRates<T> r;
    r.mu00 = mu00;
    r.mu10 = mu10;
    r.mu01 = mu01;
    r.mu11 = mu11;
    return r;
}

/**
 * The reference's fallback when the exact fit is infeasible: match E1 and E2
 * only, carrying G2 in the switching probability.
 */
template <class T>
Map<T> block_fallback(const T& E1, const T& E2, const T& G2) {
    using fitdetail::num_sqrt;
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T zero = num_traits<T>::from_int(0), four = num_traits<T>::from_int(4);
    const T SCV = T((E2 - E1 * E1) / (E1 * E1));

    if (num_traits<T>::to_double(SCV) < 1.0) {
        // SCV < 1: the reference returns the exponential; its general MAP(2)
        // branch is commented out there and is not reinstated here.
        Map<T> m;
        m.D0 = Matrix<T>(1, 1, T(-one / E1));
        m.D1 = Matrix<T>(1, 1, T(one / E1));
        return m;
    }

    // disc = 2 E1^2 (SCV - 1), so mu1 = E1 (1 - sqrt((SCV-1)/2)).
    const T disc = T(-four * E1 * E1 + two * E2);
    const T r = num_sqrt(disc > zero ? disc : zero);
    const T mu1 = T(E1 - r / two);
    const T mu2 = T(E1 + r / two);
    if (!(num_traits<T>::to_double(mu1) > 0.0))
        throw InputError(
            "map_block: the moment set is infeasible for the exact fit, and the reference's "
            "hyperexponential fallback is defined only below SCV 3, where its first branch rate "
            "E1 (1 - sqrt((SCV-1)/2)) is still positive. MATLAB returns the resulting "
            "non-generator regardless; refusing here rather than handing back a matrix that is "
            "not a MAP. Supply a feasible third moment, or fit with map_mmpp2");
    const T p = T(one / two - G2 / two);

    Matrix<T> D0(2, 2, zero);
    D0(0, 0) = T(-one / mu1);
    D0(1, 1) = T(-one / mu2);
    Matrix<T> P(2, 2, zero);
    P(0, 0) = T(one - p);
    P(0, 1) = p;
    P(1, 0) = p;
    P(1, 1) = T(one - p);

    Map<T> m;
    m.D0 = D0;
    m.D1 = Matrix<T>(2, 2, zero);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            T acc = zero;
            for (std::size_t q = 0; q < 2; ++q) acc += -D0(i, q) * P(q, j);
            m.D1(i, j) = acc;
        }
    return m;
}

}  // namespace blockdetail

/**
 * @param E1 first moment
 * @param E2 second moment
 * @param E3 third moment
 * @param G2 autocorrelation decay rate rho(i)/rho(i-1)
 */
template <class T>
Map<T> map_block(const T& E1, const T& E2, const T& E3, const T& G2) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_block inverts the moment equations through a radical");
    typedef fitdetail::Cplx<T> C;
    const T zero = num_traits<T>::from_int(0);

    const blockdetail::BlockRates<T> r =
        blockdetail::block_rates<T>(C(E1), C(E2), C(E3), C(G2));

    // The reference's own test: a residual imaginary part above 1e-4 anywhere
    // means the fit is not real and is discarded.
    const C* all[4] = {&r.mu00, &r.mu10, &r.mu01, &r.mu11};
    for (int i = 0; i < 4; ++i) {
        const double im = num_traits<T>::to_double(all[i]->im);
        if (!(std::fabs(im) <= 1e-4) || std::isnan(im)) return blockdetail::block_fallback(E1, E2, G2);
    }
    const T mu00 = r.mu00.re, mu10 = r.mu10.re, mu01 = r.mu01.re, mu11 = r.mu11.re;
    for (int i = 0; i < 4; ++i)
        if (std::isnan(num_traits<T>::to_double(all[i]->re)))
            return blockdetail::block_fallback(E1, E2, G2);
    if (mu00 < zero || mu11 < zero || mu01 < zero || mu10 < zero)
        return blockdetail::block_fallback(E1, E2, G2);

    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = T(-mu00 - mu01);
    m.D0(1, 1) = T(-mu10 - mu11);
    m.D1(0, 0) = mu00;
    m.D1(0, 1) = mu01;
    m.D1(1, 0) = mu10;
    m.D1(1, 1) = mu11;
    if (!map_isfeasible(m)) return blockdetail::block_fallback(E1, E2, G2);
    return m;
}

/** `map_block` with the SCV spelling of the second argument. */
template <class T>
Map<T> map_block_scv(const T& E1, const T& SCV, const T& E3, const T& G2) {
    return map_block(E1, T((num_traits<T>::from_int(1) + SCV) * E1 * E1), E3, G2);
}

/**
 * `map_feasblock`: repair the moments into the feasible region, then fit.
 *
 * An E2 at or below the exponential value makes the SCV non-positive, and an E3
 * below (3/2) E2^2 / E1 is outside what any MAP(2) admits; both are raised to
 * their limits plus a tolerance, exactly as the reference does.
 */
template <class T>
Map<T> map_feasblock(const T& E1, const T& E2_in, const T& E3_in, const T& G2) {
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T tol = num_traits<T>::from_double(1e-10);
    T E2 = E2_in, E3 = E3_in;

    // The exponential boundary E2 == 2 E1^2 is returned as the Poisson process,
    // rescaled to the requested mean.
    if (E2 == two * E1 * E1) {
        Map<T> m;
        m.D0 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
        m.D1 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
        m.D0(0, 0) = -one;
        m.D0(1, 1) = -one;
        m.D1(0, 0) = num_traits<T>::from_rational(1, 2);
        m.D1(0, 1) = num_traits<T>::from_rational(1, 2);
        m.D1(1, 0) = num_traits<T>::from_rational(1, 2);
        m.D1(1, 1) = num_traits<T>::from_rational(1, 2);
        return map_scale(m, E1);
    }
    if (E2 <= two * E1 * E1) E2 = T((two + tol) * E1 * E1);
    const T e3min = T(num_traits<T>::from_rational(3, 2) * E2 * E2 / E1);
    if (E3 <= e3min) E3 = T((num_traits<T>::from_rational(3, 2) + tol) * E2 * E2 / E1);
    return map_block(E1, E2, E3, G2);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_BLOCK_H
