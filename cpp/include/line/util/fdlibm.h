/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_FDLIBM_H
#define LINE_UTIL_FDLIBM_H

/**
 * The fdlibm elementary functions Java specifies, reproduced, for the code
 * paths whose whole purpose is to land on the same bits as a Java reference.
 *
 * WHY THIS FILE EXISTS: libm IS NOT A FIXED FUNCTION. `std::log1p` is allowed a
 * 1-ulp error and glibc changed which representative it returns between 2.35
 * and 2.39 (Ubuntu 22.04 vs 24.04). Measured over 200k draws, `log1p` and
 * `expm1` differ between those two while `log`, `exp`, `pow`, `sqrt`, `lgamma`
 * and `tgamma` agree. In ordinary numerics a last-bit difference is noise, but
 * a SEEDED discrete-event simulation is a chaotic map of its variates: one ulp
 * on one interarrival reorders the event queue and the whole sample path parts
 * company. That is not hypothetical -- the same `common/ldes` binary, the same
 * model.json and the same `-s 100000 --seed 23000` gave QLen 98.565592 on a
 * glibc-2.35 host and 98.562490 under the containerized MATLAB's glibc 2.39,
 * against an ABSOLUTE 1e-3 gate on the recorded baseline.
 *
 * WHY fdlibm IS THE RIGHT TARGET AND NOT MERELY A STABLE ONE. `StrictMath` IS
 * fdlibm by specification, and `Math.log1p` was measured to agree with it on
 * every one of 400k draws. Hashing raw bit patterns so the comparison is
 * language-neutral, this implementation is bit-identical to both:
 *
 *   C, host glibc 2.35        16873585727775930126
 *   C, container glibc 2.39   14249616856320735164
 *   C, this file              12005405864943986042
 *   Java Math.log1p           12005405864943986042
 *   Java StrictMath.log1p     12005405864943986042
 *
 * So routing the SSJ variate layer here does two things at once: it makes the
 * engine reproduce itself on any glibc, and it moves it ONTO the Java engine's
 * arithmetic rather than beside it. Only the seeded sample-path code should
 * call these -- the analytical APIs are free to use libm, where a last bit does
 * not cascade.
 *
 * Source: Sun's fdlibm s_log1p.c, the algorithm `StrictMath.log1p` is defined
 * to use. Transcribed with its constants and branch structure intact.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace line {
namespace fdlibm {

namespace detail {

inline int hi_word(double x) {
    std::int64_t i;
    std::memcpy(&i, &x, sizeof i);
    return static_cast<int>(i >> 32);
}

inline void set_hi_word(double* x, int hi) {
    std::int64_t i;
    std::memcpy(&i, x, sizeof i);
    i = (static_cast<std::int64_t>(hi) << 32) | (i & 0xffffffffLL);
    std::memcpy(x, &i, sizeof i);
}

}  // namespace detail

/**
 * log(1+x), fdlibm's representative -- the one StrictMath.log1p returns.
 *
 * The argument reduction writes 1+x = 2^k (1+f) with f in [sqrt(2)/2, sqrt(2)),
 * carries the rounding of that sum in the correction term c, and evaluates
 * log(1+f) from the odd series in s = f/(2+f). The magic constants are the
 * high words of the branch points: 0x3FDA827A is sqrt(2)/2 - 1 and 0x6a09e is
 * the mantissa of sqrt(2).
 */
inline double log1p(double x) {
    static const double ln2_hi = 6.93147180369123816490e-01;
    static const double ln2_lo = 1.90821492927058770002e-10;
    static const double two54 = 1.80143985094819840000e+16;
    static const double Lp1 = 6.666666666666735130e-01;
    static const double Lp2 = 3.999999999940941908e-01;
    static const double Lp3 = 2.857142874366239149e-01;
    static const double Lp4 = 2.222219843214978396e-01;
    static const double Lp5 = 1.818357216161805012e-01;
    static const double Lp6 = 1.531383769920937332e-01;
    static const double Lp7 = 1.479819860511658591e-01;

    double hfsq, f = 0.0, c = 0.0, s, z, R, u;
    int k, hx, hu = 0, ax;

    hx = detail::hi_word(x);
    ax = hx & 0x7fffffff;
    k = 1;
    if (hx < 0x3FDA827A) {          // x < 0.41422
        if (ax >= 0x3ff00000) {     // x <= -1.0
            if (x == -1.0) return -two54 / 0.0;   // -inf
            return (x - x) / (x - x);             // NaN
        }
        if (ax < 0x3e200000) {      // |x| < 2**-29
            if (two54 + x > 0.0 && ax < 0x3c900000) return x;
            return x - x * x * 0.5;
        }
        if (hx > 0 || hx <= static_cast<int>(0xbfd2bec3)) {
            k = 0;
            f = x;
            hu = 1;
        }
    }
    if (hx >= 0x7ff00000) return x + x;   // inf or NaN
    if (k != 0) {
        if (hx < 0x43400000) {
            u = 1.0 + x;
            hu = detail::hi_word(u);
            k = (hu >> 20) - 1023;
            // The correction term recovers what 1+x rounded away.
            c = (k > 0) ? 1.0 - (u - x) : x - (u - 1.0);
            c /= u;
        } else {
            u = x;
            hu = detail::hi_word(u);
            k = (hu >> 20) - 1023;
            c = 0.0;
        }
        hu &= 0x000fffff;
        if (hu < 0x6a09e) {
            detail::set_hi_word(&u, hu | 0x3ff00000);
        } else {
            k += 1;
            detail::set_hi_word(&u, hu | 0x3fe00000);
            hu = (0x00100000 - hu) >> 2;
        }
        f = u - 1.0;
    }

    hfsq = 0.5 * f * f;
    if (hu == 0) {   // |f| < 2**-20
        if (f == 0.0) {
            if (k == 0) return 0.0;
            c += k * ln2_lo;
            return k * ln2_hi + c;
        }
        R = hfsq * (1.0 - 0.66666666666666666 * f);
        if (k == 0) return f - R;
        return k * ln2_hi - ((R - (k * ln2_lo + c)) - f);
    }
    s = f / (2.0 + f);
    z = s * s;
    R = z * (Lp1 + z * (Lp2 + z * (Lp3 + z * (Lp4 + z * (Lp5 + z * (Lp6 + z * Lp7))))));
    if (k == 0) return f - (hfsq - s * (hfsq + R));
    return k * ln2_hi - ((hfsq - (s * (hfsq + R) + (k * ln2_lo + c))) - f);
}

}  // namespace fdlibm
}  // namespace line

#endif  // LINE_UTIL_FDLIBM_H
