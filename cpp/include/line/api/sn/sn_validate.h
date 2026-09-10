/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_VALIDATE_H
#define LINE_API_SN_SN_VALIDATE_H

/**
 * Consistency checks on a NetworkStruct, a port of
 * jar/src/main/java/jline/api/sn/SnValidate.java and ValidationLevel.java.
 *
 * The reference returns a LIST OF MESSAGES rather than throwing, so a caller can
 * report every defect of a struct at once instead of the first one; that
 * contract is kept here. The setters in sn_setters.h are the intended callers:
 * they mutate a struct in place and can leave it internally inconsistent, and
 * this is what tells them so.
 *
 * WHAT DIFFERS FROM THE JAR, and why. The reference tests the SHAPE of the
 * matrices it stores (rates, scv, nservers, njobs, classprio, phases each being
 * an independently-sized Matrix that could disagree with nstations/nclasses).
 * This port holds nservers, the populations, the priorities and the phase counts
 * as per-station and per-class members of Station and JobClass, so those shapes
 * cannot disagree by construction and the corresponding checks would be vacuous.
 * The checks that survive are the ones with real content: the two genuinely
 * free-standing matrices (rates, scv), the parallel `disabled` mask that
 * replaces MATLAB's NaN sentinel, and the value checks on rates, populations,
 * routing and servers.
 *
 * NaN IS NOT A DISABLED MARKER HERE. The reference skips a NaN rate because
 * MATLAB writes NaN for a pair the class never visits; this port carries that
 * out of band in `disabled` (Rational has no NaN), so the skip reads the flag.
 *
 * ARITHMETIC: field plus comparisons, so it instantiates under Rational.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/** How much of the struct to check; the reference's ValidationLevel. */
enum class ValidationLevel {
    Full,     ///< every check
    Minimal,  ///< dimensions only
    None      ///< skip everything
};

namespace detail {

inline std::string idx2(std::size_t i, std::size_t j) {
    return "[" + std::to_string(i) + "," + std::to_string(j) + "]";
}

}  // namespace detail

/** The (nstations x nclasses) matrices agree with the declared dimensions. */
template <class T>
std::vector<std::string> sn_validate_dimensions(const qn::NetworkStruct<T>& sn) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::vector<std::string> errors;
    const std::string want =
        " do not match (nstations=" + std::to_string(M) + " x nclasses=" + std::to_string(K) + ")";
    if (sn.rates.rows() != M || sn.rates.cols() != K)
        errors.push_back("rates matrix dimensions (" + std::to_string(sn.rates.rows()) + "x" +
                         std::to_string(sn.rates.cols()) + ")" + want);
    if (sn.scv.rows() != M || sn.scv.cols() != K)
        errors.push_back("scv matrix dimensions (" + std::to_string(sn.scv.rows()) + "x" +
                         std::to_string(sn.scv.cols()) + ")" + want);
    if (!sn.disabled.empty()) {
        if (sn.disabled.size() != M)
            errors.push_back("disabled has " + std::to_string(sn.disabled.size()) +
                             " rows, expected nstations=" + std::to_string(M));
        for (std::size_t i = 0; i < sn.disabled.size(); ++i)
            if (sn.disabled[i].size() != K) {
                errors.push_back("disabled row " + std::to_string(i) + " has " +
                                 std::to_string(sn.disabled[i].size()) +
                                 " entries, expected nclasses=" + std::to_string(K));
                break;
            }
    }
    if (sn.stations.size() != M)
        errors.push_back("stations list holds " + std::to_string(sn.stations.size()) +
                         " entries, expected nstations=" + std::to_string(M));
    if (sn.classes.size() != K)
        errors.push_back("classes list holds " + std::to_string(sn.classes.size()) +
                         " entries, expected nclasses=" + std::to_string(K));
    return errors;
}

/** Rates and SCVs are non-negative wherever the pair is enabled. */
template <class T>
std::vector<std::string> sn_validate_rates(const qn::NetworkStruct<T>& sn) {
    std::vector<std::string> errors;
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < sn.rates.rows(); ++i)
        for (std::size_t j = 0; j < sn.rates.cols(); ++j) {
            if (i < sn.disabled.size() && j < sn.disabled[i].size() && sn.disabled[i][j]) continue;
            if (sn.rates(i, j) < zero)
                errors.push_back("rates" + detail::idx2(i, j) + " = " +
                                 std::to_string(num_traits<T>::to_double(sn.rates(i, j))) +
                                 " is negative");
        }
    for (std::size_t i = 0; i < sn.scv.rows(); ++i)
        for (std::size_t j = 0; j < sn.scv.cols(); ++j) {
            if (i < sn.disabled.size() && j < sn.disabled[i].size() && sn.disabled[i][j]) continue;
            if (sn.scv(i, j) < zero)
                errors.push_back("scv" + detail::idx2(i, j) + " = " +
                                 std::to_string(num_traits<T>::to_double(sn.scv(i, j))) +
                                 " is negative (SCV must be >= 0)");
        }
    return errors;
}

/** Class populations are neither NaN nor negative; an open class may be Inf. */
template <class T>
std::vector<std::string> sn_validate_population(const qn::NetworkStruct<T>& sn) {
    std::vector<std::string> errors;
    for (std::size_t k = 0; k < sn.classes.size(); ++k) {
        const double n = sn.classes[k].population;
        if (std::isnan(n))
            errors.push_back("njobs for class " + std::to_string(k) + " is NaN");
        else if (std::isfinite(n) && n < 0)
            errors.push_back("njobs for class " + std::to_string(k) + " = " + std::to_string(n) +
                             " is negative");
    }
    return errors;
}

/** Every non-empty row of rt is a probability vector summing to one. */
template <class T>
std::vector<std::string> sn_validate_routing(const qn::NetworkStruct<T>& sn) {
    std::vector<std::string> errors;
    const T zero = num_traits<T>::from_int(0);
    const double tolerance = 1e-6;
    for (std::size_t i = 0; i < sn.rt.rows(); ++i) {
        T rowSum = zero;
        bool hasNonZero = false;
        for (std::size_t j = 0; j < sn.rt.cols(); ++j) {
            const T v = sn.rt(i, j);
            if (v < zero)
                errors.push_back("rt" + detail::idx2(i, j) + " = " +
                                 std::to_string(num_traits<T>::to_double(v)) + " is negative");
            if (v > zero) hasNonZero = true;
            rowSum += v;
        }
        const double s = num_traits<T>::to_double(rowSum);
        if (hasNonZero && std::fabs(s - 1.0) > tolerance)
            errors.push_back("rt row " + std::to_string(i) + " sum = " + std::to_string(s) +
                             " (expected 1.0)");
    }
    return errors;
}

/** Server counts are positive, except at a Source or a Sink. */
template <class T>
std::vector<std::string> sn_validate_servers(const qn::NetworkStruct<T>& sn) {
    std::vector<std::string> errors;
    for (std::size_t i = 0; i < sn.stations.size(); ++i) {
        const double n = sn.stations[i].nservers;
        if (std::isnan(n)) {
            errors.push_back("nservers[" + std::to_string(i) + "] is NaN");
        } else if (n <= 0 && !std::isinf(n)) {
            const qn::NodeType nt = sn.stations[i].nodetype;
            if (nt != qn::NodeType::Source && nt != qn::NodeType::Sink)
                errors.push_back("nservers[" + std::to_string(i) + "] = " + std::to_string(n) +
                                 " must be positive");
        }
    }
    return errors;
}

/**
 * @param sn    the struct to check
 * @param level how much of it to check
 * @return one message per defect, empty when the struct is consistent
 */
template <class T>
std::vector<std::string> sn_validate(const qn::NetworkStruct<T>& sn,
                                     ValidationLevel level = ValidationLevel::Full) {
    std::vector<std::string> errors;
    if (level == ValidationLevel::None) return errors;
    const std::vector<std::string> dims = sn_validate_dimensions(sn);
    errors.insert(errors.end(), dims.begin(), dims.end());
    if (level != ValidationLevel::Full) return errors;
    const std::vector<std::string> r = sn_validate_rates(sn);
    errors.insert(errors.end(), r.begin(), r.end());
    const std::vector<std::string> p = sn_validate_population(sn);
    errors.insert(errors.end(), p.begin(), p.end());
    const std::vector<std::string> rt = sn_validate_routing(sn);
    errors.insert(errors.end(), rt.begin(), rt.end());
    const std::vector<std::string> ns = sn_validate_servers(sn);
    errors.insert(errors.end(), ns.begin(), ns.end());
    return errors;
}

// ---- index guards, returned as a message or the empty string ------------

template <class T>
std::string sn_validate_station_index(const qn::NetworkStruct<T>& sn, std::size_t ist,
                                      const std::string& paramName = "stationIdx") {
    if (ist >= sn.nstations)
        return paramName + "=" + std::to_string(ist) + " is out of bounds [0, " +
               std::to_string(sn.nstations == 0 ? 0 : sn.nstations - 1) + "]";
    return std::string();
}

template <class T>
std::string sn_validate_class_index(const qn::NetworkStruct<T>& sn, std::size_t r,
                                    const std::string& paramName = "classIdx") {
    if (r >= sn.nclasses)
        return paramName + "=" + std::to_string(r) + " is out of bounds [0, " +
               std::to_string(sn.nclasses == 0 ? 0 : sn.nclasses - 1) + "]";
    return std::string();
}

template <class T>
std::string sn_validate_node_index(const qn::NetworkStruct<T>& sn, std::size_t ind,
                                   const std::string& paramName = "nodeIdx") {
    const std::size_t n = sn.nodes.size();
    if (ind >= n)
        return paramName + "=" + std::to_string(ind) + " is out of bounds [0, " +
               std::to_string(n == 0 ? 0 : n - 1) + "]";
    return std::string();
}

template <class T>
std::string sn_validate_node_type(const qn::NetworkStruct<T>& sn, std::size_t ind,
                                  qn::NodeType expected) {
    if (ind >= sn.nodes.size())
        return "nodeIdx=" + std::to_string(ind) + " is out of bounds";
    const qn::NodeType actual = sn.nodes[ind].nodetype;
    if (actual != expected)
        return "Node " + std::to_string(ind) + " is " + lang::node_type_to_text(actual) +
               ", expected " + lang::node_type_to_text(expected);
    return std::string();
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_VALIDATE_H
