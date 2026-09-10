/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_NUM_NUMBER_H
#define LINE_NUM_NUMBER_H

/**
 * Number-type abstraction for the templated API port.
 *
 * Three arithmetic modes are exposed to callers:
 *   double  - IEEE 754, what MATLAB / the JAR / native Python use today
 *   exact   - arbitrary-precision rationals, no rounding at all
 *   real    - fixed high-precision binary floating point
 *
 * The default backends are header-only Boost.Multiprecision types
 * (Boost Software License 1.0), so the default build carries no LGPL
 * obligation and a redistributable Python wheel remains possible. Defining
 * LINE_MP_USE_GMP switches the exact and real backends to GMP mpq and MPFR,
 * which are faster but LGPL; the algorithms are unchanged, only the typedef.
 *
 * Every algorithm is a template on the number type T and consults
 * num_traits<T> for what T can do. Algorithms that need log/exp/pow assert
 * num_traits<T>::has_transcendental at compile time, so instantiating an
 * inherently inexact algorithm at exact arithmetic is a build error rather
 * than a silent fallback.
 */

#include <cmath>
#include <string>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#ifdef LINE_MP_USE_GMP
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>
#endif

namespace line {

// Backend typedefs

/* `Rational` IS PINNED TO et_off ON PURPOSE, and the `et_off` is the whole
 * point of spelling these out rather than using Boost's `cpp_rational` /
 * `mpq_rational` aliases, which are `et_on`.
 *
 * With expression templates on, `a + b` does not evaluate: it builds an
 * `expression<...>` node holding REFERENCES to its operands. That is a
 * use-after-free the moment the node outlives them, and a DEDUCED return type
 * is how it escapes -- `[](Rational x) { return Rational(4) * x - Rational(1); }`
 * returns a node over two temporaries that die with the return statement.
 * It cost two SIGSEGVs (test_pfqn_mvac_oi_clw.cpp 2026-08-25,
 * test_rootfind.cpp 2026-08-26), and both were invisible on the 22.04
 * development host: Boost 1.72 added rvalue-reference handling to the
 * operators, so 1.74 collapses that expression eagerly while 1.71 (every 20.04
 * worker) does not. doctest cannot catch a SIGSEGV, so each one aborted the
 * binary mid-run and left several hundred cases with no verdict at all.
 *
 * Turning the templates off removes the entire class rather than the two sites
 * that happened to be found. It cannot change a single digit: rational
 * arithmetic here is exact, so expression templates were only ever eliding
 * temporaries, never altering a value. Do NOT "restore" the Boost aliases.
 *
 * BigInt and Real are left alone: `Real` is already et_off for these backends,
 * and no BigInt expression is returned through a deduced type. */
#ifdef LINE_MP_USE_GMP
using Rational =
    boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off>;
using BigInt = boost::multiprecision::mpz_int;
template <unsigned Digits10>
using Real = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<Digits10>>;
#else
using Rational = boost::multiprecision::number<boost::multiprecision::cpp_rational_backend,
                                               boost::multiprecision::et_off>;
using BigInt = boost::multiprecision::cpp_int;
template <unsigned Digits10>
using Real = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<Digits10>>;
#endif

/** Precision tiers offered by the CLI's --arith real:`<digits>` flag. */
using Real50 = Real<50>;
using Real100 = Real<100>;
using Real200 = Real<200>;

// ---------------------------------------------------------------------------
// log of a big integer, without ever forming the value as a double
// ---------------------------------------------------------------------------

/**
 * log(v) for a positive arbitrary-precision integer. Keeps only the leading
 * 53 bits of the mantissa and accounts for the discarded bits in the exponent,
 * so the result is finite for values far outside the double range.
 */
inline double log_bigint(const BigInt& v) {
    if (v == 0) return -std::numeric_limits<double>::infinity();
    BigInt a = v < 0 ? BigInt(-v) : v;
    long bits = static_cast<long>(boost::multiprecision::msb(a)) + 1;
    long shift = bits > 53 ? bits - 53 : 0;
    BigInt mant = a >> shift;
    return std::log(static_cast<double>(mant)) + static_cast<double>(shift) * std::log(2.0);
}

// ---------------------------------------------------------------------------
// num_traits
// ---------------------------------------------------------------------------

template <class T>
struct num_traits;  // intentionally undefined for unsupported types

template <>
struct num_traits<double> {
    using type = double;
    static constexpr bool is_exact = false;
    static constexpr bool has_transcendental = true;
    static const char* name() { return "double"; }

    static double from_int(long v) { return static_cast<double>(v); }
    static double from_rational(long num, long den) {
        return static_cast<double>(num) / static_cast<double>(den);
    }
    static double from_double(double v) { return v; }
    static double to_double(const double& v) { return v; }
    /** log of the value, always returned as a double. */
    static double log_as_double(const double& v) { return std::log(v); }
    static std::string to_string(const double& v) { return std::to_string(v); }
};

template <>
struct num_traits<Rational> {
    using type = Rational;
    static constexpr bool is_exact = true;
    static constexpr bool has_transcendental = false;
    static const char* name() { return "exact"; }

    static Rational from_int(long v) { return Rational(v); }
    static Rational from_rational(long num, long den) { return Rational(num, den); }
    /** Exact: a double is a dyadic rational, so this conversion loses nothing. */
    static Rational from_double(double v) { return Rational(v); }
    static double to_double(const Rational& v) { return static_cast<double>(v); }
    static double log_as_double(const Rational& v) {
        return log_bigint(BigInt(numerator(v))) - log_bigint(BigInt(denominator(v)));
    }
    static std::string to_string(const Rational& v) { return v.str(); }
    static std::string numerator_str(const Rational& v) { return BigInt(numerator(v)).str(); }
    static std::string denominator_str(const Rational& v) { return BigInt(denominator(v)).str(); }
};

template <unsigned D>
struct num_traits<Real<D>> {
    using type = Real<D>;
    static constexpr bool is_exact = false;
    static constexpr bool has_transcendental = true;
    static constexpr unsigned digits10 = D;
    static const char* name() { return "real"; }

    static Real<D> from_int(long v) { return Real<D>(v); }
    static Real<D> from_rational(long num, long den) { return Real<D>(num) / Real<D>(den); }
    static Real<D> from_double(double v) { return Real<D>(v); }
    static double to_double(const Real<D>& v) { return static_cast<double>(v); }
    static double log_as_double(const Real<D>& v) { return static_cast<double>(log(v)); }
    static std::string to_string(const Real<D>& v) { return v.str(); }
};

// ---------------------------------------------------------------------------
// Generic helpers usable from any algorithm
// ---------------------------------------------------------------------------

template <class T>
inline T num_abs(const T& v) {
    using std::abs;
    return abs(v);
}

template <>
inline double num_abs<double>(const double& v) {
    return std::fabs(v);
}

/** Factorial as a value of T. Exact for Rational and for BigInt-backed types. */
template <class T>
inline T num_factorial(unsigned n) {
    T f = num_traits<T>::from_int(1);
    for (unsigned k = 2; k <= n; ++k) f *= num_traits<T>::from_int(static_cast<long>(k));
    return f;
}

/** Integer power, valid in any field (no transcendental requirement). */
template <class T>
inline T num_pow_int(const T& base, unsigned e) {
    T r = num_traits<T>::from_int(1);
    T b = base;
    unsigned k = e;
    while (k > 0) {
        if (k & 1u) r *= b;
        b *= b;
        k >>= 1;
    }
    return r;
}

}  // namespace line

#endif  // LINE_NUM_NUMBER_H
