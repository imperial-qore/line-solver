/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMPP2_FIT3_H
#define LINE_API_MAM_MMPP2_FIT3_H

/**
 * MMPP(2) matching three moments and the autocorrelation decay rate
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fit3.m).
 *
 * G2 is the constant ratio rho(i)/rho(i-1) of consecutive lag
 * autocorrelations. The reference is a symbolic-toolbox solution of the
 * moment-matching system, so the four parameters mu00, mu11, q01, q10 are
 * single rational expressions in E1, E3, SCV and G2 over the moment
 * discriminant. They are transcribed verbatim here, with the three repeated
 * subexpressions named:
 *
 *   DISC = the discriminant polynomial under the square root,
 *   SQ   = sqrt(DISC),
 *   pw(x, k) = x^k.
 *
 * A vanishing G2 degenerates to an uncorrelated MAP(1)-style fit (mu11 = 0),
 * which the reference handles with its own closed form; that branch is
 * reproduced exactly, including the g2tol threshold below which it fires.
 *
 * The resulting representation is
 *   D0 = [ -mu00-q01  q01 ; q10  -mu11-q10 ],  D1 = diag(mu00, mu11).
 *
 * Gated on transcendental arithmetic: SQ is a square root of a moment
 * discriminant with no exact rational counterpart.
 */

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * MMPP(2) with moments (E1, E2, E3) and autocorrelation decay rate G2.
 * g2tol is the threshold below which the uncorrelated branch is taken
 * (MATLAB uses 1e-6).
 */
template <class T>
Map<T> mmpp2_fit3(const T& E1, const T& E2, const T& E3, const T& G2, const T& g2tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fit3 requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    if (E1 == zero) throw InputError("mmpp2_fit3: zero first moment");
    const T SCV = (E2 - E1 * E1) / (E1 * E1);

    T mu00, mu11, q01, q10;
    if (G2 < g2tol || G2 == zero) {
        mu00 = 2*(6*pw(E1,3)*SCV-E3)/E1/(6*pw(E1,3)*SCV+3*pw(E1,3)*pw(SCV,2)+3*pw(E1,3)-2*E3);
        mu11 = zero;
        q01 =  9*pw(E1,5)*(SCV-1)*(pw(SCV,2)-2*SCV+1)/(6*pw(E1,3)*SCV-E3)/(6*pw(E1,3)*SCV+3*pw(E1,3)*pw(SCV,2)+
              3*pw(E1,3)-2*E3);
        q10 = -3*(SCV-1)*pw(E1,2)/(6*pw(E1,3)*SCV-E3);
    } else {
        const T DISC = (pw(E3,2)-12*pw(E1,3)*SCV*E3+6*pw(E1,3)*G2*E3-6*G2*SCV*pw(E1,3)*E3+18*G2*pw(SCV,3)*pw(E1,6)-
                          18*pw(E1,6)*G2*pw(SCV,2)+9*pw(E1,6)*pw(G2,2)+36*pw(E1,6)*pw(SCV,2)+18*pw(E1,6)*G2*SCV-18*pw(E1,6)*
                          SCV*pw(G2,2)+9*pw(E1,6)*pw(SCV,2)*pw(G2,2)-18*pw(E1,6)*G2);
        if (DISC < zero) throw NumericError("mmpp2_fit3: negative moment discriminant");
        const T SQ = num_sqrt(DISC);
        mu00 = G2*(-4*E3*G2+4*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-
               6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3*G2-18*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
               SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2-18*pw(E1,3)*(-3*pw(E1,3)*
               G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
               2*E3)*G2*pw(SCV,2)-12*pw(E1,3)*pw(G2,2)-12*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
               SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,2)*SCV+12*pw(E1,3)*
               (-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-
               3*pw(E1,3)+2*E3)*G2*SCV+12*pw(E1,3)*G2*pw(SCV,2)-9*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-
               6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*SCV+3*pw(E1,3)*
               (-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-
               3*pw(E1,3)+2*E3)+12*pw(E1,3)*pw(G2,2)*SCV+9*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
               SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(SCV,2)+12*pw(E1,3)*G2+
               12*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-
               6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,2)-3*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
               SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(SCV,3))/(12*pw(E1,3)*pw(G2,3)*
               SCV+3*pw(E1,3)*pw(SCV,3)*G2-12*pw(E1,3)*pw(G2,3)+18*pw(E1,3)*pw(G2,2)*pw(SCV,2)-3*pw(E1,3)*
               G2+27*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-
               6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2*pw(SCV,2)-9*pw(E1,3)*G2*pw(SCV,2)+18*pw(E1,3)*pw(G2,2)-12*
               pw(E1,3)*pw(G2,2)*SCV+9*pw(E1,3)*G2*SCV-12*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
               SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,3)*SCV-9*pw(E1,3)*(-
               3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-
               3*pw(E1,3)+2*E3)*pw(SCV,3)*G2-24*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+
               E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,2)*pw(SCV,2)-(-3*pw(E1,3)*
               G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
               2*E3)*E3*pw(SCV,2)+4*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-
               6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3*pw(G2,2)+12*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*
               pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,3)-(-3*pw(E1,3)*
               G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
               2*E3)*E3+2*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*
               pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3*SCV+9*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
               SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2+24*pw(E1,3)*(-3*pw(E1,3)*
               G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
               2*E3)*pw(G2,2)*SCV-27*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*
               pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2*SCV+6*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*
               G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*SCV-12*
               pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*
               SCV-3*pw(E1,3)+2*E3)*pw(SCV,2)-24*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+
               E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,2)+6*pw(E1,3)*(-3*pw(E1,3)*
               G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
               2*E3)*pw(SCV,3)-4*E3*pw(G2,2))/E1;
        mu11 = (-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/E1/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*
               SCV-3*pw(E1,3)+2*E3);
        q01 = -3*pw(E1,2)*(-6*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*
              pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*SCV+12*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
              SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2*SCV-6*G2*SCV*
              pw(E1,2)-3*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-
              6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2+(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/
              E1/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3+3*pw(E1,2)*G2+6*(-3*pw(E1,3)*G2+
              3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
              2*E3)*pw(SCV,2)-9*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*
              pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(SCV,2)*G2+3*pw(E1,2)*G2*pw(SCV,2)-E3*(-3*pw(E1,3)*
              G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/E1/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
              2*E3)*SCV-6*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-
              6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,2)*SCV+6*pw(E1,2)*pw(G2,2)*SCV+3*(-3*pw(E1,3)*G2+3*pw(E1,3)*
              G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*
              pw(G2,2)-G2*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/E1/(-3*pw(E1,3)*pw(SCV,2)-
              6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3-3*pw(E1,2)*pw(G2,2)+3*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-
              6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(SCV,2)*
              pw(G2,2)-3*pw(E1,2)*pw(SCV,2)*pw(G2,2)+G2*SCV*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
              SCV+E3+SQ)/E1/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3)/(-45*(-3*pw(E1,3)*
              G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*
              pw(E1,3)+2*E3)*G2*pw(SCV,2)+18*pw(G2,2)*pw(E1,5)*SCV+18*pw(E1,5)*pw(G2,3)-27*pw(E1,5)*pw(G2,2)*
              pw(SCV,2)+6*pw(E1,2)*pw(G2,2)*E3-27*pw(E1,5)*pw(G2,2)-18*pw(E1,5)*pw(G2,3)*SCV-18*pw(E1,5)*
              G2*SCV+18*pw(E1,5)*G2*pw(SCV,2)+3*pw(E1,2)*G2*E3-3*pw(E1,2)*G2*E3*SCV+(-3*pw(E1,3)*G2+3*pw(E1,3)*
              G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/E1/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(E3,2)+
              3*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-6*
              pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2*SCV*E3-36*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+
              E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,2)*SCV+36*(-3*
              pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*
              SCV-3*pw(E1,3)+2*E3)*pw(G2,2)+36*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/
              (-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(SCV,2)+45*(-3*pw(E1,3)*G2+3*pw(E1,3)*
              G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*
              G2*SCV-12*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-
              6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*SCV*E3-3*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+
              SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2*E3+9*(-3*pw(E1,3)*G2+
              3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+
              2*E3)*G2*pw(SCV,3)+36*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*
              pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,2)*pw(SCV,2)-6*(-3*pw(E1,3)*G2+3*pw(E1,3)*
              G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,2)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*
              pw(G2,2)*E3+18*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*
              pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,3)*SCV-18*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-
              6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(G2,3)-
              9*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)*pw(E1,5)/(-3*pw(E1,3)*pw(SCV,2)-6*
              pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2);
        q10 = 3*(-3*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-
              6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(SCV,3)-3*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*
              SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2*pw(SCV,2)+6*pw(E1,3)*pw(SCV,2)+
              3*pw(E1,3)*G2*pw(SCV,2)+3*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/
              (-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*pw(SCV,2)+6*pw(E1,3)*(-3*pw(E1,3)*G2+
              3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*
              E3)*G2*SCV-E3*SCV-6*pw(E1,3)*SCV+(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-
              3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3*SCV-6*pw(E1,3)*G2*SCV-3*pw(E1,3)*(-
              3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-
              3*pw(E1,3)+2*E3)*SCV-(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-
              6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*E3+3*pw(E1,3)*G2-3*pw(E1,3)*(-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-
              6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-3*pw(E1,3)+2*E3)*G2+3*pw(E1,3)*
              (-3*pw(E1,3)*G2+3*pw(E1,3)*G2*SCV-6*pw(E1,3)*SCV+E3+SQ)/(-3*pw(E1,3)*pw(SCV,2)-6*pw(E1,3)*SCV-
              3*pw(E1,3)+2*E3)+E3)*pw(E1,2)*(-1+G2)/DISC;
    }

    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = -mu00 - q01;
    m.D0(0, 1) = q01;
    m.D0(1, 0) = q10;
    m.D0(1, 1) = -mu11 - q10;
    m.D1(0, 0) = mu00;
    m.D1(1, 1) = mu11;
    return m;
}

/** mmpp2_fit3 with the MATLAB default g2tol = 1e-6. */
template <class T>
Map<T> mmpp2_fit3(const T& E1, const T& E2, const T& E3, const T& G2) {
    return mmpp2_fit3(E1, E2, E3, G2, T(num_traits<T>::from_double(1e-6)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMPP2_FIT3_H
