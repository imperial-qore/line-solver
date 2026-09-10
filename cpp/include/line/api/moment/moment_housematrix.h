/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_HOUSEMATRIX_H
#define LINE_API_MOMENT_MOMENT_HOUSEMATRIX_H

/**
 * Conversion matrix of one edge of the house of moments.
 *
 * Templated port of matlab/src/api/moment/moment_housematrix.m. The house of
 * moments collects the raw, factorial, upper-factorial, binomial,
 * negative-binomial and tail sequences of a discrete law; every edge of it is a
 * linear map with an integer or rational table, so the exact instantiation
 * returns the transform with no rounding at all.
 *
 * The tail edges are the only UPPER triangular ones, so a caller must apply the
 * table with a full matrix-vector product and not with the lower-triangular
 * apply_table helper.
 *
 * Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 */

#include <string>
#include <vector>

#include "line/api/moment/moment_lah.h"
#include "line/api/moment/moment_stirling1.h"
#include "line/api/moment/moment_stirling2.h"
#include "line/api/moment/moment_stirlingcycle.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace moment {

/** Edge labels of the house of moments, one per MATLAB edge string. */
enum class MomentEdge {
    FactorialFromRaw,
    RawFromFactorial,
    UpfactorialFromRaw,
    RawFromUpfactorial,
    BinomialFromFactorial,
    NegbinomialFromUpfactorial,
    FactorialFromBinomial,
    UpfactorialFromNegbinomial,
    FactorialFromUpfactorial,
    UpfactorialFromFactorial,
    NegbinomialFromBinomial,
    BinomialFromNegbinomial,
    BinomialFromTail,
    TailFromBinomial
};

/** Parses the MATLAB edge string into its label. */
inline MomentEdge moment_edge_from_string(const std::string& edge) {
    if (edge == "factorial_from_raw") return MomentEdge::FactorialFromRaw;
    if (edge == "raw_from_factorial") return MomentEdge::RawFromFactorial;
    if (edge == "upfactorial_from_raw") return MomentEdge::UpfactorialFromRaw;
    if (edge == "raw_from_upfactorial") return MomentEdge::RawFromUpfactorial;
    if (edge == "binomial_from_factorial") return MomentEdge::BinomialFromFactorial;
    if (edge == "negbinomial_from_upfactorial") return MomentEdge::NegbinomialFromUpfactorial;
    if (edge == "factorial_from_binomial") return MomentEdge::FactorialFromBinomial;
    if (edge == "upfactorial_from_negbinomial") return MomentEdge::UpfactorialFromNegbinomial;
    if (edge == "factorial_from_upfactorial") return MomentEdge::FactorialFromUpfactorial;
    if (edge == "upfactorial_from_factorial") return MomentEdge::UpfactorialFromFactorial;
    if (edge == "negbinomial_from_binomial") return MomentEdge::NegbinomialFromBinomial;
    if (edge == "binomial_from_negbinomial") return MomentEdge::BinomialFromNegbinomial;
    if (edge == "binomial_from_tail") return MomentEdge::BinomialFromTail;
    if (edge == "tail_from_binomial") return MomentEdge::TailFromBinomial;
    throw InputError("moment_housematrix: unknown edge " + edge);
}

/** (n+1)x(n+1) conversion table of the given edge. */
template <class T>
Matrix<T> moment_housematrix(MomentEdge edge, int n) {
    if (n < 0)
        throw InputError("moment_housematrix: the maximum order n must be a nonnegative integer");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    Matrix<T> Tm(static_cast<std::size_t>(n) + 1, static_cast<std::size_t>(n) + 1, zero);
    switch (edge) {
        case MomentEdge::FactorialFromRaw:
            return moment_stirling1<T>(n);
        case MomentEdge::RawFromFactorial:
            return moment_stirling2<T>(n);
        case MomentEdge::UpfactorialFromRaw:
            return moment_stirlingcycle<T>(n);
        case MomentEdge::RawFromUpfactorial: {
            Matrix<T> S = moment_stirling2<T>(n);
            for (int i = 0; i <= n; ++i)
                for (int j = 0; j <= i; ++j)
                    Tm(i, j) = ((i - j) % 2 == 0) ? S(i, j) : -S(i, j);
            return Tm;
        }
        case MomentEdge::BinomialFromFactorial:
        case MomentEdge::NegbinomialFromUpfactorial:
            for (int i = 0; i <= n; ++i)
                Tm(i, i) = one / num_factorial<T>(static_cast<unsigned>(i));
            return Tm;
        case MomentEdge::FactorialFromBinomial:
        case MomentEdge::UpfactorialFromNegbinomial:
            for (int i = 0; i <= n; ++i) Tm(i, i) = num_factorial<T>(static_cast<unsigned>(i));
            return Tm;
        case MomentEdge::FactorialFromUpfactorial:
        case MomentEdge::UpfactorialFromFactorial: {
            Matrix<T> L = moment_lah<T>(n);
            Tm(0, 0) = one;
            for (int i = 1; i <= n; ++i)
                for (int k = 1; k <= i; ++k)
                    Tm(i, k) = (edge == MomentEdge::UpfactorialFromFactorial || (i - k) % 2 == 0)
                                   ? L(i, k)
                                   : -L(i, k);
            return Tm;
        }
        case MomentEdge::NegbinomialFromBinomial:
        case MomentEdge::BinomialFromNegbinomial:
            Tm(0, 0) = one;
            for (int i = 1; i <= n; ++i)
                for (int k = 1; k <= i; ++k) {
                    const T c = num_nck<T>(i - 1, k - 1);
                    Tm(i, k) = (edge == MomentEdge::NegbinomialFromBinomial || (i - k) % 2 == 0)
                                   ? c
                                   : -c;
                }
            return Tm;
        case MomentEdge::BinomialFromTail:
        case MomentEdge::TailFromBinomial:
            Tm(0, 0) = one;
            for (int i = 1; i <= n; ++i)
                for (int k = i; k <= n; ++k) {
                    const T c = num_nck<T>(k - 1, i - 1);
                    Tm(i, k) =
                        (edge == MomentEdge::BinomialFromTail || (k - i) % 2 == 0) ? c : -c;
                }
            return Tm;
    }
    throw InputError("moment_housematrix: unknown edge");
}

/** String overload matching the MATLAB call signature. */
template <class T>
Matrix<T> moment_housematrix(const std::string& edge, int n) {
    return moment_housematrix<T>(moment_edge_from_string(edge), n);
}

/** Full matrix-vector product, needed because the tail edges are upper triangular. */
template <class T>
std::vector<T> moment_apply_full(const Matrix<T>& A, const std::vector<T>& v) {
    std::vector<T> r(v.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < v.size(); ++j) r[i] += A(i, j) * v[j];
    return r;
}

}  // namespace moment
}  // namespace line

#endif
