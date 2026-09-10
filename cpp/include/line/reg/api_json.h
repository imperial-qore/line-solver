/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_REG_API_JSON_H
#define LINE_REG_API_JSON_H

/**
 * The one conversion policy between JSON and the templated API layer.
 *
 * Every host boundary that carries API arguments and results as JSON goes
 * through this header: the --api path of the CLI today, a pybind11 or MEX
 * gateway later. Two gateways that each invent their own number conversion
 * will disagree the first time a caller writes 0.6, so the policy is stated
 * once, here, and both sides share it.
 *
 * ARGUMENTS. The object's keys are the MATLAB parameter names verbatim
 * ({"L": [[0.6,0.4]], "N": [2,1], "Z": [1,0.5]}). A 2-D array is a matrix,
 * row-major with the outer index the row; a 1-D array is a row vector; a bare
 * number is a 1x1 scalar. An unrecognised key is an error, never ignored: a
 * misspelt argument that silently takes its default is a wrong answer.
 *
 * NUMBERS. A JSON number is interpreted as the SHORTEST DECIMAL LITERAL that
 * round-trips to it, and that decimal is what the arithmetic sees. So 0.6
 * becomes the rational 3/5 in exact arithmetic, not the dyadic
 * 5404319552844595/9007199254740992 that double-to-rational conversion would
 * give. This is the only reading under which "exact" means what a caller
 * writing 0.6 intends, and it is deterministic because the shortest
 * round-tripping decimal of a double is unique. A caller who wants a value
 * that has no short decimal form passes a STRING: "1/3" and "0.3333" are both
 * accepted, and the string is taken literally at every arithmetic.
 *
 * RESULTS. A value of the algorithm's number type T encodes as
 *   double     -> the bare JSON number
 *   exact      -> {"double": `<approx>`, "num": "`<decimal>`", "den": "`<decimal>`"}
 *   real:`<D>`   -> {"double": `<approx>`, "dec": "`<decimal string>`"}
 * Exact numerators and denominators overflow every integer type, so they cross
 * as decimal strings; the host rebuilds them with sym(num)/sym(den) or
 * fractions.Fraction. A value that is a C++ double in the algorithm itself
 * regardless of T -- lG is the only one in the port -- always encodes as a
 * bare JSON number, because there is no exact value to report: lG is computed
 * exponent-safely as log(num) - log(den) and is finite where G is not
 * representable at all. Never exponentiate it back.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>

#include "json.hpp"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace reg {

using Json = nlohmann::json;

// ---------------------------------------------------------------------------
// Decimal literal recovery and exact decimal parsing
// ---------------------------------------------------------------------------

/**
 * The shortest decimal literal that round-trips to v. nlohmann::json does not
 * retain the text of a number, so this reconstructs the literal a human wrote:
 * for every double that came from a short decimal in the file, the reconstruction
 * IS that decimal, because a double has at most one shortest round-tripping
 * representation.
 */
inline std::string shortest_decimal(double v) {
    char buf[64];
    for (int prec = 1; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof(buf), "%.*g", prec, v);
        if (std::strtod(buf, nullptr) == v) return std::string(buf);
    }
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

inline BigInt pow10_bigint(unsigned e) {
    BigInt p = 1;
    for (unsigned k = 0; k < e; ++k) p *= 10;
    return p;
}

/**
 * Exact value of a decimal literal, with no rounding anywhere: sign, digits,
 * optional fraction and optional exponent are read symbolically and assembled
 * as numerator over a power of ten.
 */
inline Rational rational_from_decimal(const std::string& text) {
    const std::size_t slash = text.find('/');
    if (slash != std::string::npos)
        return rational_from_decimal(text.substr(0, slash)) /
               rational_from_decimal(text.substr(slash + 1));

    std::size_t i = 0;
    const std::size_t n = text.size();
    while (i < n && (text[i] == ' ' || text[i] == '\t')) ++i;
    bool negative = false;
    if (i < n && (text[i] == '+' || text[i] == '-')) negative = (text[i++] == '-');

    std::string digits;
    long frac_digits = 0;
    bool any = false;
    while (i < n && text[i] >= '0' && text[i] <= '9') {
        digits += text[i++];
        any = true;
    }
    if (i < n && text[i] == '.') {
        ++i;
        while (i < n && text[i] >= '0' && text[i] <= '9') {
            digits += text[i++];
            ++frac_digits;
            any = true;
        }
    }
    if (!any) throw InputError("not a decimal number: '" + text + "'");

    long exponent = 0;
    if (i < n && (text[i] == 'e' || text[i] == 'E')) {
        ++i;
        bool eneg = false;
        if (i < n && (text[i] == '+' || text[i] == '-')) eneg = (text[i++] == '-');
        std::string edig;
        while (i < n && text[i] >= '0' && text[i] <= '9') edig += text[i++];
        if (edig.empty()) throw InputError("malformed exponent in '" + text + "'");
        exponent = std::strtol(edig.c_str(), nullptr, 10);
        if (eneg) exponent = -exponent;
    }
    while (i < n && (text[i] == ' ' || text[i] == '\t')) ++i;
    if (i != n) throw InputError("trailing characters in number '" + text + "'");

    // leading-zero-as-octal fix: see _kb/14-cpp-multiprecision.md
    std::size_t first = digits.find_first_not_of('0');
    digits = (first == std::string::npos) ? std::string("0") : digits.substr(first);
    const BigInt mantissa(digits);
    Rational r(mantissa);
    const long scale = exponent - frac_digits;
    if (scale > 0)
        r *= Rational(pow10_bigint(static_cast<unsigned>(scale)));
    else if (scale < 0)
        r /= Rational(pow10_bigint(static_cast<unsigned>(-scale)));
    return negative ? Rational(-r) : r;
}

/** double value of a decimal literal, including the "a/b" fraction form. */
inline double double_from_decimal(const std::string& text) {
    const std::size_t slash = text.find('/');
    if (slash != std::string::npos)
        return double_from_decimal(text.substr(0, slash)) /
               double_from_decimal(text.substr(slash + 1));
    const char* s = text.c_str();
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    if (end == s) throw InputError("not a decimal number: '" + text + "'");
    while (*end == ' ' || *end == '\t') ++end;
    if (*end != '\0') throw InputError("trailing characters in number '" + text + "'");
    return v;
}

// ---------------------------------------------------------------------------
// Decimal literal -> the algorithm's number type
// ---------------------------------------------------------------------------

template <class T>
struct NumFromDecimal;

template <>
struct NumFromDecimal<double> {
    static double parse(const std::string& text) { return double_from_decimal(text); }
};

template <>
struct NumFromDecimal<Rational> {
    static Rational parse(const std::string& text) { return rational_from_decimal(text); }
};

template <unsigned D>
struct NumFromDecimal<Real<D>> {
    static Real<D> parse(const std::string& text) {
        const std::size_t slash = text.find('/');
        if (slash != std::string::npos)
            return Real<D>(parse(text.substr(0, slash))) / Real<D>(parse(text.substr(slash + 1)));
        // Reject what the backend would accept loosely, so a typo is an error
        // rather than a silent zero, then let the backend do the rounding.
        (void)double_from_decimal(text);
        return Real<D>(text);
    }
};

/** The decimal literal behind a JSON scalar, per the policy in the file header. */
inline std::string decimal_text(const Json& j, const std::string& where) {
    if (j.is_string()) return j.get<std::string>();
    if (j.is_number_integer()) return std::to_string(j.get<long long>());
    if (j.is_number_unsigned()) return std::to_string(j.get<unsigned long long>());
    if (j.is_number_float()) return shortest_decimal(j.get<double>());
    if (j.is_boolean()) return j.get<bool>() ? "1" : "0";
    throw InputError(where + ": expected a number or a numeric string, got " +
                     std::string(j.type_name()));
}

template <class T>
T number_from_json(const Json& j, const std::string& where) {
    return NumFromDecimal<T>::parse(decimal_text(j, where));
}

// ---------------------------------------------------------------------------
// Shapes
// ---------------------------------------------------------------------------

/**
 * A matrix from a JSON value: 2-D array as rows, 1-D array as a row vector,
 * scalar as 1x1, empty array as the empty matrix. Ragged rows are an error.
 */
template <class T>
Matrix<T> matrix_from_json(const Json& j, const std::string& where) {
    if (!j.is_array()) {
        Matrix<T> m(1, 1);
        m(0, 0) = number_from_json<T>(j, where);
        return m;
    }
    if (j.empty()) return Matrix<T>();
    if (j[0].is_array()) {
        const std::size_t rows = j.size();
        const std::size_t cols = j[0].size();
        for (std::size_t i = 0; i < rows; ++i) {
            if (!j[i].is_array())
                throw InputError(where + ": row " + std::to_string(i) + " is not an array");
            if (j[i].size() != cols)
                throw InputError(where + ": row " + std::to_string(i) + " has " +
                                 std::to_string(j[i].size()) + " entries, row 0 has " +
                                 std::to_string(cols));
        }
        if (cols == 0) return Matrix<T>();
        Matrix<T> m(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t k = 0; k < cols; ++k)
                m(i, k) = number_from_json<T>(j[i][k], where);
        return m;
    }
    Matrix<T> m(1, j.size());
    for (std::size_t k = 0; k < j.size(); ++k) m(0, k) = number_from_json<T>(j[k], where);
    return m;
}

/** A flat numeric vector: 1-D array, or a 1-row / 1-column 2-D array. */
template <class T>
std::vector<T> vector_from_json(const Json& j, const std::string& where) {
    const Matrix<T> m = matrix_from_json<T>(j, where);
    if (m.empty()) return std::vector<T>();
    if (m.rows() != 1 && m.cols() != 1)
        throw InputError(where + ": expected a vector, got a " + std::to_string(m.rows()) + "x" +
                         std::to_string(m.cols()) + " matrix");
    std::vector<T> v(m.size());
    for (std::size_t k = 0; k < m.size(); ++k) v[k] = m[k];
    return v;
}

/** An integer vector; a non-integral entry is an error, never a truncation. */
inline std::vector<int> int_vector_from_json(const Json& j, const std::string& where) {
    std::vector<Json> flat;
    if (!j.is_array()) {
        flat.push_back(j);
    } else {
        for (const Json& e : j) {
            if (e.is_array())
                for (const Json& f : e) flat.push_back(f);
            else
                flat.push_back(e);
        }
    }
    std::vector<int> v;
    v.reserve(flat.size());
    for (const Json& e : flat) {
        const std::string text = decimal_text(e, where);
        const Rational r = rational_from_decimal(text);
        if (denominator(r) != 1)
            throw InputError(where + ": expected an integer, got " + text);
        const double d = static_cast<double>(r);
        if (d > 2147483647.0 || d < -2147483648.0)
            throw InputError(where + ": integer out of range: " + text);
        v.push_back(static_cast<int>(d));
    }
    return v;
}

// ---------------------------------------------------------------------------
// Encoding results
// ---------------------------------------------------------------------------

template <class T>
struct EncodeScalar;

template <>
struct EncodeScalar<double> {
    static Json encode(const double& v) { return Json(v); }
};

template <>
struct EncodeScalar<Rational> {
    static Json encode(const Rational& v) {
        Json j;
        j["double"] = num_traits<Rational>::to_double(v);
        j["num"] = num_traits<Rational>::numerator_str(v);
        j["den"] = num_traits<Rational>::denominator_str(v);
        return j;
    }
};

template <unsigned D>
struct EncodeScalar<Real<D>> {
    static Json encode(const Real<D>& v) {
        Json j;
        j["double"] = num_traits<Real<D>>::to_double(v);
        j["dec"] = v.str(static_cast<std::streamsize>(D), std::ios_base::scientific);
        return j;
    }
};

template <class T>
Json encode_scalar(const T& v) {
    return EncodeScalar<T>::encode(v);
}

template <class T>
Json encode_vector(const std::vector<T>& v) {
    Json a = Json::array();
    for (const T& x : v) a.push_back(encode_scalar(x));
    return a;
}

template <class T>
Json encode_matrix(const Matrix<T>& m) {
    Json a = Json::array();
    for (std::size_t i = 0; i < m.rows(); ++i) {
        Json row = Json::array();
        for (std::size_t k = 0; k < m.cols(); ++k) row.push_back(encode_scalar(m(i, k)));
        a.push_back(row);
    }
    return a;
}

inline Json encode_ints(const std::vector<int>& v) {
    Json a = Json::array();
    for (int x : v) a.push_back(x);
    return a;
}

/** A list of matrices, the shape the MMAP/BMAP families return as {D0,D1,...}. */
template <class T>
Json encode_matrices(const std::vector<Matrix<T> >& v) {
    Json a = Json::array();
    for (const Matrix<T>& m : v) a.push_back(encode_matrix(m));
    return a;
}

/** A list of vectors, the shape a per-class or per-segment result returns. */
template <class T>
Json encode_vectors(const std::vector<std::vector<T> >& v) {
    Json a = Json::array();
    for (const std::vector<T>& x : v) a.push_back(encode_vector(x));
    return a;
}

/**
 * A count. Written as a plain JSON integer at every arithmetic: an iteration
 * count or a dimension is exact in all of them, so wrapping it in the
 * {"double",...} envelope the scalars use would suggest a precision question
 * that a cardinal number does not have.
 */
inline Json encode_count(std::size_t n) { return Json(static_cast<std::uint64_t>(n)); }

// ---------------------------------------------------------------------------
// Argument object
// ---------------------------------------------------------------------------

/**
 * Named-argument reader over the parsed JSON object. Every read marks the key;
 * done() then refuses any key the function did not ask for, naming it and
 * listing what the function does accept. Silently ignoring an unknown key is
 * how a caller ends up with the default value of the argument they thought
 * they were setting.
 */
class Args {
public:
    Args(const Json& j, std::string function) : j_(j), fn_(std::move(function)) {
        if (!j_.is_object())
            throw InputError(fn_ + ": the argument JSON must be an object keyed by the MATLAB "
                                   "parameter names, got " +
                             std::string(j_.type_name()));
    }

    bool has(const char* key) const { return j_.find(key) != j_.end(); }

    /** Required argument; throws when absent. */
    const Json& get(const char* key) {
        seen_.insert(key);
        auto it = j_.find(key);
        if (it == j_.end()) throw InputError(fn_ + ": missing required argument '" + key + "'");
        return *it;
    }

    /** Optional argument; returns nullptr when absent. */
    const Json* opt(const char* key) {
        seen_.insert(key);
        auto it = j_.find(key);
        return it == j_.end() ? nullptr : &*it;
    }

    template <class T>
    Matrix<T> matrix(const char* key) {
        return matrix_from_json<T>(get(key), fn_ + ": " + key);
    }

    template <class T>
    Matrix<T> matrix_or_empty(const char* key) {
        const Json* v = opt(key);
        return v ? matrix_from_json<T>(*v, fn_ + ": " + key) : Matrix<T>();
    }

    template <class T>
    std::vector<T> vector_or_empty(const char* key) {
        const Json* v = opt(key);
        return v ? vector_from_json<T>(*v, fn_ + ": " + key) : std::vector<T>();
    }

    /** Required vector argument; absent is an error, not an empty vector. */
    template <class T>
    std::vector<T> vector(const char* key) {
        return vector_from_json<T>(get(key), fn_ + ": " + key);
    }

    /**
     * Required list of matrices, e.g. the per-class D1c of an MMAP. A single
     * matrix is NOT silently promoted to a one-element list: the two shapes
     * mean different models and the caller must say which one it meant.
     */
    template <class T>
    std::vector<Matrix<T> > matrices(const char* key) {
        const Json& j = get(key);
        const std::string where = fn_ + ": " + key;
        if (!j.is_array())
            throw InputError(where + ": expected a list of matrices, got " +
                             std::string(j.type_name()));
        std::vector<Matrix<T> > out;
        for (std::size_t i = 0; i < j.size(); ++i)
            out.push_back(matrix_from_json<T>(j[i], where + "[" + std::to_string(i) + "]"));
        return out;
    }

    template <class T>
    std::vector<Matrix<T> > matrices_or_empty(const char* key) {
        seen_.insert(key);
        if (!has(key)) return std::vector<Matrix<T> >();
        return matrices<T>(key);
    }

    /** Required numeric argument, at the arithmetic of the call. */
    template <class T>
    T number(const char* key) {
        return number_from_json<T>(get(key), fn_ + ": " + key);
    }

    /**
     * A required nonnegative count. Refuses a negative value naming the
     * argument rather than wrapping it around into a huge unsigned.
     */
    std::size_t count(const char* key) {
        const int v = required_integer(key);
        if (v < 0)
            throw InputError(fn_ + ": '" + key + "' is a count and cannot be negative, got " +
                             std::to_string(v));
        return static_cast<std::size_t>(v);
    }

    std::size_t count(const char* key, std::size_t fallback) {
        if (!has(key)) {
            seen_.insert(key);
            return fallback;
        }
        return count(key);
    }

    /** The same as count(), for the arguments the port declares `unsigned`. */
    unsigned uinteger(const char* key) { return static_cast<unsigned>(count(key)); }

    unsigned uinteger(const char* key, unsigned fallback) {
        return static_cast<unsigned>(count(key, static_cast<std::size_t>(fallback)));
    }

    int required_integer(const char* key) {
        const std::vector<int> got = int_vector_from_json(get(key), fn_ + ": " + key);
        if (got.size() != 1) throw InputError(fn_ + ": '" + key + "' must be a single integer");
        return got[0];
    }

    /** A flag. Accepts a JSON boolean or the integers 0 and 1, nothing else. */
    bool boolean(const char* key, bool fallback) {
        const Json* v = opt(key);
        if (!v) return fallback;
        if (v->is_boolean()) return v->get<bool>();
        const std::vector<int> got = int_vector_from_json(*v, fn_ + ": " + key);
        if (got.size() != 1 || (got[0] != 0 && got[0] != 1))
            throw InputError(fn_ + ": '" + key + "' must be true, false, 0 or 1");
        return got[0] != 0;
    }

    /** A string-valued option, e.g. a method or convention name. */
    std::string text(const char* key, const std::string& fallback) {
        const Json* v = opt(key);
        if (!v) return fallback;
        if (!v->is_string())
            throw InputError(fn_ + ": '" + key + "' must be a string, got " +
                             std::string(v->type_name()));
        return v->get<std::string>();
    }

    std::vector<int> ints(const char* key) {
        return int_vector_from_json(get(key), fn_ + ": " + key);
    }

    std::vector<int> ints_or_empty(const char* key) {
        const Json* v = opt(key);
        return v ? int_vector_from_json(*v, fn_ + ": " + key) : std::vector<int>();
    }

    int integer(const char* key, int fallback) {
        const Json* v = opt(key);
        if (!v) return fallback;
        const std::vector<int> got = int_vector_from_json(*v, fn_ + ": " + key);
        if (got.size() != 1)
            throw InputError(fn_ + ": '" + key + "' must be a single integer");
        return got[0];
    }

    template <class T>
    T scalar(const char* key, const T& fallback) {
        const Json* v = opt(key);
        return v ? number_from_json<T>(*v, fn_ + ": " + key) : fallback;
    }

    /**
     * Refuse an argument the JSON boundary cannot carry, naming it, rather than
     * proceeding as if it had not been given.
     */
    void unsupported(const char* key, const std::string& why) {
        seen_.insert(key);
        if (j_.find(key) != j_.end())
            throw UnsupportedError(fn_ + ": argument '" + key + "' cannot cross the JSON API "
                                                                "boundary: " +
                                   why);
    }

    void done() const {
        for (auto it = j_.begin(); it != j_.end(); ++it) {
            if (seen_.find(it.key()) != seen_.end()) continue;
            std::string accepted;
            for (const std::string& k : seen_) {
                if (!accepted.empty()) accepted += ", ";
                accepted += k;
            }
            throw InputError(fn_ + ": unknown argument '" + it.key() + "'; " + fn_ +
                             " accepts: " + accepted);
        }
    }

private:
    const Json& j_;
    std::string fn_;
    std::set<std::string> seen_;
};

}  // namespace reg
}  // namespace line

#endif  // LINE_REG_API_JSON_H
