/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_MOMENT_SIMPLE_H
#define LINE_API_TRACE_MTRACE_MOMENT_SIMPLE_H

/**
 * Class-pair cross moments of a marked trace.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_moment_simple.m,
 * cross-checked against
 * jar/src/main/java/jline/api/trace/Mtrace_moment_simple.java.
 *
 * In BOTH codebases this function is a byte-for-byte duplicate of
 * mtrace_cross_moment (same body, same doc comment, only the name differs),
 * so it is kept as an alias rather than a second implementation: any future
 * change to one of them would otherwise silently diverge here.
 *
 * ARITHMETIC: as mtrace_cross_moment, exact in Rational.
 */

#include <vector>

#include "line/api/trace/mtrace_cross_moment.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"

namespace line {
namespace trace {

/** @see mtrace_cross_moment */
template <class T>
MtraceCrossMomentResult<T> mtrace_moment_simple(const std::vector<T>& Tv,
                                                const std::vector<int>& L, unsigned k) {
    return mtrace_cross_moment(Tv, L, k);
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_MOMENT_SIMPLE_H
