/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_NUM_COMPLEX_NUMBER_H
#define LINE_NUM_COMPLEX_NUMBER_H

/**
 * `std::complex<double>` as a number type for the generic linear algebra.
 *
 * WHY THIS EXISTS. One algorithm in the port needs complex arithmetic: the
 * Laplace-domain transient QBD (`mam_transient2`, `mam_transient2_open`),
 * whose level blocks are shifted by `-s I` at a COMPLEX quadrature node s and
 * whose fundamental matrices G and R are therefore complex. Everything it needs
 * -- products, identity, inverse, integer powers, LU with partial pivoting --
 * already exists in `util/linalg.h` and `util/lu.h` templated on the element
 * type, so declaring the traits here makes the whole stack available at complex
 * argument instead of duplicating it.
 *
 * WHAT IS DELIBERATELY NARROW. `to_double` returns the REAL PART, because the
 * only consumer is a Laplace inversion whose quadrature sums `Re(eta_k F(s_k))`
 * and discards the imaginary part by construction. `log_as_double` is the log
 * of the MODULUS. Neither is a faithful "convert to double"; both are the
 * meaning the transient QBD needs, and no other algorithm instantiates this
 * type. Do not reach for it as a general complex facility without revisiting
 * those two.
 *
 * `is_exact` is false and `has_transcendental` is true, so any algorithm gated
 * on exact arithmetic refuses at complex argument, as it should.
 */

#include <cmath>
#include <complex>
#include <string>

#include "line/num/number.h"
#include "line/util/lu.h"

namespace line {

using Complex = std::complex<double>;

template <>
struct num_traits<Complex> {
    using type = Complex;
    static constexpr bool is_exact = false;
    static constexpr bool has_transcendental = true;
    static const char* name() { return "complex"; }

    static Complex from_int(long v) { return Complex(static_cast<double>(v), 0.0); }
    static Complex from_rational(long num, long den) {
        return Complex(static_cast<double>(num) / static_cast<double>(den), 0.0);
    }
    static Complex from_double(double v) { return Complex(v, 0.0); }
    /** The real part: the Laplace quadrature keeps Re and discards Im. */
    static double to_double(const Complex& v) { return v.real(); }
    /** Log of the modulus. */
    static double log_as_double(const Complex& v) { return std::log(std::abs(v)); }
    static std::string to_string(const Complex& v) {
        return std::to_string(v.real()) + (v.imag() < 0.0 ? "-" : "+") +
               std::to_string(std::fabs(v.imag())) + "i";
    }
};

/** Partial pivoting compares moduli, since the complex field is unordered. */
template <class R>
struct pivot_mag<std::complex<R>> {
    using type = R;
    static R of(const std::complex<R>& x) { return std::abs(x); }
};

}  // namespace line

#endif  // LINE_NUM_COMPLEX_NUMBER_H
