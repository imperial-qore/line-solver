/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_ILTCME_TABLE_H
#define LINE_API_MAM_ILTCME_TABLE_H

/**
 * The vendored concentrated-matrix-exponential (CME) coefficient table that
 * `matlab_ilt` reads from `iltcme.json`.
 *
 * IT IS A .cpp, NOT A HEADER, AND THAT IS THE POINT. The table is ~1.2 MB of
 * double literals. As an `inline constexpr` array in a header it would be
 * re-parsed by every translation unit that includes it, and `line_cli.cpp`
 * already includes every solver header; compiling it once into `line_mp_api`
 * costs one TU instead. This header carries only the declarations, so the port
 * stays header-only for everything except this one data blob.
 *
 * VENDORED THIRD-PARTY DATA whose licence terms are NOT STATED upstream. See
 * THIRD-PARTY-NOTICES.md, which records the decision to vendor with that known.
 */

#include <cstddef>

namespace line {
namespace mam {
namespace iltcme {

/**
 * One CME entry, carrying only the fields `matlab_ilt` reads.
 *
 * `iltcme.json` also has `optim`, `phi`, `lognorm` and `mu2`; none is used by
 * the inverse transform, so none is vendored.
 */
struct CmeEntry {
    int n;              ///< number of cosine/sine terms; the transform costs n+1 evaluations
    double c;           ///< constant term
    double omega;       ///< angular frequency
    double mu1;         ///< first moment scale
    double cv2;         ///< squared coefficient of variation; smaller is steeper
    const double* a;    ///< cosine coefficients, length n
    const double* b;    ///< sine coefficients, length n
    std::size_t len;    ///< n, as a size for bounds checks
};

/** The reachable entries, in table order. */
extern const CmeEntry kTable[];
extern const std::size_t kTableSize;

}  // namespace iltcme
}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_ILTCME_TABLE_H
