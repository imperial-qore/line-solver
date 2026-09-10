/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_DECIMAL_H
#define LINE_UTIL_DECIMAL_H

/**
 * Decimal literal -> T, without a detour through double when T is exact.
 *
 * This is the boundary at which a model file becomes numbers, and it is the
 * one place where "exact arithmetic" has to decide what it is exact ABOUT. A
 * host demand written `0.01` in an .lqnx file denotes the rational 1/100. The
 * double nearest to it is 0.01000000000000000020816681711721685...; routing
 * that literal through `num_traits<Rational>::from_double` would make the
 * exact backend compute, with no rounding at all, the answer to a model that
 * is not the one on disk. So the exact backend parses the decimal digits
 * directly into num/10^k, and the inexact backends keep strtod, which is what
 * MATLAB's str2double and Java's Double.parseDouble also do.
 *
 * The consequence to keep in mind when comparing the two runs: they are not
 * two evaluations of one arithmetic problem, they are the exact solution of
 * the declared model versus the floating-point solution of its double
 * rounding. Their difference is bounded below by the input rounding, of order
 * 1e-17 relative here, and any larger gap is accumulated error in the solver.
 */

#include <cctype>
#include <cstdlib>
#include <string>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {

namespace detail {

/** Split a decimal literal into (sign, digits, exponent) with value = sign * digits * 10^exp. */
inline bool decimal_parts(const std::string& s, bool& neg, std::string& digits, long& exp10) {
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    neg = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) neg = (s[i++] == '-');

    digits.clear();
    exp10 = 0;
    bool any = false;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
        digits.push_back(s[i++]);
        any = true;
    }
    if (i < s.size() && s[i] == '.') {
        ++i;
        while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
            digits.push_back(s[i++]);
            --exp10;
            any = true;
        }
    }
    if (!any) return false;
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        ++i;
        bool eneg = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) eneg = (s[i++] == '-');
        long e = 0;
        bool anye = false;
        while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
            e = e * 10 + (s[i++] - '0');
            anye = true;
        }
        if (!anye) return false;
        exp10 += eneg ? -e : e;
    }
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    return i == s.size();
}

}  // namespace detail

/** The literal as an exact rational, num/10^k with no rounding. */
inline Rational rational_from_decimal(const std::string& s) {
    bool neg = false;
    std::string digits;
    long exp10 = 0;
    if (!detail::decimal_parts(s, neg, digits, exp10))
        throw InputError("num_from_decimal: '" + s + "' is not a decimal literal");
    BigInt mant = 0;
    for (char c : digits) mant = mant * 10 + (c - '0');
    BigInt num = mant;
    BigInt den = 1;
    if (exp10 >= 0) {
        for (long k = 0; k < exp10; ++k) num *= 10;
    } else {
        for (long k = 0; k < -exp10; ++k) den *= 10;
    }
    Rational r(num, den);
    return neg ? Rational(-r) : r;
}

/**
 * Parse a decimal literal into T.
 *
 * Rational reconstructs num/10^k from the digits; every other backend uses
 * strtod, matching what MATLAB's str2double and Java's Double.parseDouble do.
 */
template <class T>
T num_from_decimal(const std::string& s) {
    return num_traits<T>::from_double(std::strtod(s.c_str(), nullptr));
}

template <>
inline Rational num_from_decimal<Rational>(const std::string& s) {
    return rational_from_decimal(s);
}

/** Parse a decimal literal as a plain double (multiplicities, populations, tolerances). */
inline double dbl_from_decimal(const std::string& s, double fallback) {
    if (s.empty()) return fallback;
    const char* p = s.c_str();
    char* end = nullptr;
    const double v = std::strtod(p, &end);
    if (end == p) return fallback;
    return v;
}

}  // namespace line

#endif  // LINE_UTIL_DECIMAL_H
