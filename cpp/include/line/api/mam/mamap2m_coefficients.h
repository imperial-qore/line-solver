/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAMAP2M_COEFFICIENTS_H
#define LINE_API_MAM_MAMAP2M_COEFFICIENTS_H

/**
 * The marking coefficients of a canonical AMAP(2), for the sigma fitters.
 *
 * Templated port of matlab/lib/m3a/m3a/mamap2m/mamap2m_can1_coefficients.m and
 * mamap2m_can2_coefficients.m, cross-checked against
 * jar/src/main/java/jline/api/mam/Mamap2m_coefficients.java.
 *
 * `mamap22_fit_fs_multiclass` and `mamap22_fit_bs_multiclass` match the class
 * TRANSITION probabilities (sigma) alongside a forward or backward moment, and
 * the relation between the marking and those characteristics is not the simple
 * affine one the F+B fitter uses. These two tables carry it: G (or E) holds the
 * per-flow contributions to the class probability, the sigma and the moment; U
 * (or V) the quadratic terms; and Y (or Z) three determinants of G that decide
 * whether the system is solvable.
 *
 * REFERENCE TYPO, corrected here and already corrected in the JAR. MATLAB's
 * `mamap2m_can1_coefficients.m` assigns `G(10)` TWICE in consecutive lines:
 *
 *     G(10) = (r1*r2^2)/(r1*r2 - r2 + 1);      % this is G(9)
 *     G(10) = h1 - (h1*r1)/(r2*(r1 - 1) + 1);
 *
 * so the first value is discarded and G(9) is left at zero. The JAR writes them
 * to indices 9 and 10 respectively, which is the only reading under which the
 * table has no hole, and that is what is done here. Y does not read G(9), so the
 * three determinants are unaffected; a consumer that reads the ninth
 * coefficient gets zero from MATLAB and the intended value here.
 *
 * ARITHMETIC: field. Rational expressions only, so this instantiates exactly.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** The three coefficient tables of one canonical form. */
template <class T>
struct Mamap2mCoefficients {
    std::vector<T> G;  ///< 15 entries for form 1, 14 for form 2 (there called E)
    std::vector<T> U;  ///< 12 entries (V for form 2)
    std::vector<T> Y;  ///< 3 determinants (Z for form 2)
};

/**
 * First canonical form, a positive autocorrelation decay.
 *
 * Indices are 1-based in the reference; the vectors here are 0-based, so
 * G[0] is the reference's G(1).
 */
template <class T>
Mamap2mCoefficients<T> mamap2m_can1_coefficients(const T& h1, const T& h2, const T& r1,
                                                 const T& r2) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2), three = num_traits<T>::from_int(3);
    Mamap2mCoefficients<T> c;
    c.G.assign(15, zero);
    c.U.assign(12, zero);
    c.Y.assign(3, zero);

    const T A = T(r2 * (r1 - one) + one);   // the reference's r2*(r1-1)+1
    const T B = T(r1 * r2 - r2 + one);      // and r1*r2 - r2 + 1
    if (A == zero || B == zero)
        throw NumericError(
            "mamap2m_can1_coefficients: the canonical denominators vanish at this (r1, r2); the "
            "marking coefficients are undefined there");

    c.G[0] = T(one - r1 / A);
    c.G[1] = T(-(r1 * (r2 - one)) / B);
    c.G[2] = T((r1 * r2) / B);
    c.G[3] = T((r1 * (r1 - one)) / A - r1 + one);
    c.G[4] = T(-(r1 * (r1 - one) * (r2 - one) * (r2 - two)) / B);
    c.G[5] = T((r1 * r2 * (r1 - one) * (r2 - one)) / B);
    c.G[6] = T((r1 * r1 * (r2 - one) * (r2 - one)) / A);
    c.G[7] = T(-(r1 * r2 * (r1 + one) * (r2 - one)) / B);
    // The reference writes this into G(10) and immediately overwrites it; see
    // the header note. It belongs at G(9).
    c.G[8] = T((r1 * r2 * r2) / B);
    c.G[9] = T(h1 - (h1 * r1) / A);
    c.G[10] = T(-(r1 * (r2 - one) * (h1 + h2 - h1 * r2)) / B);
    c.G[11] = T((r1 * r2 * (h1 + h2 - h1 * r2)) / B);
    c.G[12] = T(((h1 + h2 * r1) * (r1 - one) * (r2 - one)) / B);
    c.G[13] = T(-(r1 * (h1 + h2 * r1) * (r2 - one)) / B);
    c.G[14] = T((h2 * r1 * r2) / B);

    const T t1 = T(h1 - h2 + h2 * r1);
    c.U[0] = T(B * B);
    c.U[1] = T(-B * (two * h1 - h1 * r1 - two * h1 * r2 + three * h2 * r1 - h2 * r1 * r1 +
                     h2 * r1 * r1 * r2 + h1 * r1 * r2 - h2 * r1 * r2));
    c.U[2] = T(r1 * (r2 - one) * t1 * t1);
    c.U[3] = T(B * (h2 * h2 * r1 - h1 * h1 * r2 + h1 * h1 + h1 * h2 * r1 - h1 * h2 * r1 * r2));
    c.U[4] = T(-r1 * (r2 - one) * B * t1);
    c.U[5] = T(r1 * (r2 - one) * (h1 - h1 * r2 + h2 * r1) * t1);
    const T t2 = T(h2 - h1 * r2);
    c.U[6] = T(B * B);
    c.U[7] = T(-B * (two * h1 - two * h1 * r2 + h2 * r1 - h1 * r1 * r2 * r2 + h1 * r1 * r2 +
                     h2 * r1 * r2));
    c.U[8] = T(r1 * t2 * t2 * (r2 - one));
    c.U[9] = T(B * (h2 * h2 * r1 - h1 * h1 * r2 + h1 * h1 + h1 * h2 * r1 - h1 * h2 * r1 * r2));
    c.U[10] = T(-r1 * t2 * (r2 - one) * B);
    c.U[11] = T(r1 * t2 * (r2 - one) * (h1 - h1 * r2 + h2 * r1));

    // Y(1) is the 3x3 determinant of the (G1,G2,G3 | G10,G11,G12 | G13,G14,G15)
    // block, and Y(2), Y(3) two of its 2x2 minors.
    c.Y[0] = T(c.G[0] * c.G[10] * c.G[14] - c.G[0] * c.G[11] * c.G[13] -
               c.G[1] * c.G[9] * c.G[14] + c.G[1] * c.G[11] * c.G[12] +
               c.G[2] * c.G[9] * c.G[13] - c.G[2] * c.G[10] * c.G[12]);
    c.Y[1] = T(c.G[2] * c.G[12] - c.G[0] * c.G[14]);
    c.Y[2] = T(c.G[9] * c.G[2] - c.G[11] * c.G[0]);
    return c;
}

/** Second canonical form, a negative autocorrelation decay (E, V, Z). */
template <class T>
Mamap2mCoefficients<T> mamap2m_can2_coefficients(const T& h1, const T& h2, const T& r1,
                                                 const T& r2) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2), three = num_traits<T>::from_int(3);
    Mamap2mCoefficients<T> c;
    c.G.assign(14, zero);
    c.U.assign(12, zero);
    c.Y.assign(3, zero);

    const T A = T(r2 * (r1 - one) - r1 + two);
    const T B = T(r1 * (r2 - one) - r2 + two);
    const T C = T(r1 + r2 - r1 * r2 - two);
    if (A == zero || B == zero || C == zero)
        throw NumericError(
            "mamap2m_can2_coefficients: the canonical denominators vanish at this (r1, r2); the "
            "marking coefficients are undefined there");

    c.G[0] = T(one - one / A);
    c.G[1] = T(-(r2 - one) / B);
    c.G[2] = T(r2 / B);
    c.G[3] = T((r2 - two) / B - r2 + two);
    c.G[4] = T(r2 - r2 / B);
    c.G[5] = T(-(r1 * (r2 - one) * (r2 - one)) / C);
    c.G[6] = T(-r2 - (r2 * (two * r2 - three)) / B);
    c.G[7] = T(r2 * r2 / A);
    c.G[8] = T(h1 - h1 / A);
    c.G[9] = T(h1 * (r2 - one) - ((r2 - one) * (two * h1 + h2 - h1 * r2)) / B);
    c.G[10] = T((r2 * (two * h1 + h2 - h1 * r2)) / B - h1 * r2);
    c.G[11] = T(h2 - h2 / A);
    c.G[12] = T(((h1 + h2 * r1) * (r2 - one)) / C);
    c.G[13] = T((h2 * r2) / B);

    const T t1 = T(h1 - h2 + h2 * r1);
    const T t2 = T(h1 - h2 - h1 * r1 + h1 * r1 * r2);
    c.U[0] = T(-C * C);
    c.U[1] = T(-C * (two * h1 + two * h2 - h1 * r2 - h2 * r2 + h2 * r1 * r2));
    c.U[2] = T(h2 * (two * h1 - h1 * r2 + h2 * r1) * C);
    c.U[3] = T((r2 - one) * t1 * t1);
    c.U[4] = T(t1 * (two * r2 - r1 * r2 + r1 * r2 * r2 - r2 * r2));
    c.U[5] = T(-(h1 * r2 + h2 * r2 - h1 * r2 * r2) * t1);
    c.U[6] = T(-C * C);
    c.U[7] = T(-C * (two * h1 + two * h2 - h1 * r2 - h2 * r2 + h1 * r1 * r2 * r2 - h1 * r1 * r2));
    c.U[8] = T(h1 * C * (two * h2 + h1 * r1 - h2 * r2 + h1 * r1 * r2 * r2 - two * h1 * r1 * r2));
    c.U[9] = T((r2 - one) * t2 * t2);
    c.U[10] = T(-r2 * t2 * C);
    c.U[11] = T(-r2 * (h1 + h2 - h1 * r2) * t2);

    c.Y[0] = T(c.G[9] * c.G[11] * c.G[2] - c.G[9] * c.G[13] * c.G[0] -
               c.G[10] * c.G[11] * c.G[1] + c.G[10] * c.G[12] * c.G[0] -
               c.G[12] * c.G[2] * c.G[8] + c.G[13] * c.G[1] * c.G[8]);
    c.Y[1] = T(c.G[11] * c.G[1] - c.G[12] * c.G[0]);
    c.Y[2] = T(c.G[9] * c.G[0] - c.G[1] * c.G[8]);
    return c;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAMAP2M_COEFFICIENTS_H
