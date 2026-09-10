/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_RNG_SSJ_H
#define LINE_UTIL_RNG_SSJ_H

/**
 * The two random number generators the Java LDES engine draws from, reproduced
 * exactly: SSJ's `MRG32k3a` and `java.util.Random`.
 *
 * WHY THIS EXISTS. `common/ldes` is to become the C++ engine compiled, in place
 * of a GraalVM image of the Java one, and a simulator that answers with a
 * different sample path at the same seed is not a replacement -- every seeded
 * golden in MATLAB, Python and the parity harness would re-baseline, and
 * "LDES agrees across the codebases at a given seed", which `CLAUDE.md` states
 * as an invariant, would weaken to "agrees in distribution". The C++ engine
 * previously drew from `std::mt19937_64`; this header is the first half of
 * closing that gap, and `ldes_sampler.h` is the second (the variate algorithms
 * must consume the same uniforms in the same order).
 *
 * WHAT MADE IT TRACTABLE. `Solver_ssj` never relies on SSJ's package-level
 * stream sequence: every stream is created and then given an EXPLICIT six-long
 * seed derived from the run seed and a per-(node, class) offset, e.g.
 *
 *     long offset = ((long) (numSources + svcIdx) * numClasses + k) * 10 + 1000;
 *     stream.setSeed(new long[] { seed + offset, ..., seed + offset + 5 });
 *
 * so none of SSJ's 2^127 jump-ahead machinery is on the path. Matching the
 * engine means matching the recurrence, `nextValue`/`nextDouble`, and those
 * offsets -- not the stream hierarchy.
 *
 * MRG32k3a is L'Ecuyer's combined multiple recursive generator (Operations
 * Research 47(1), 1999): two order-three recurrences modulo m1 = 2^32 - 209 and
 * m2 = 2^32 - 22853, combined by subtraction. SSJ computes it in DOUBLES with
 * the multipliers split so every intermediate stays exactly representable, and
 * that arithmetic is reproduced here rather than an integer rewrite: the double
 * form is what the reference runs, and the two agree only if the rounding does.
 *
 * `java.util.Random` is the specified 48-bit LCG of the Java Language
 * Specification: seed scrambling by 0x5DEECE66D, `next(bits)` taking the high
 * bits, `nextDouble` from a 26-bit and a 27-bit draw, and `nextInt(bound)` with
 * its rejection loop for non-power-of-two bounds. All of it is normative, so a
 * faithful transcription is exact by construction.
 */

#include <cmath>
#include <cstdint>
#include <string>

#include "line/util/error.h"

namespace line {
namespace rng {

// ---------------------------------------------------------------------------
// SSJ MRG32k3a
// ---------------------------------------------------------------------------

/**
 * SSJ's `umontreal.ssj.rng.MRG32k3a`, state and all.
 *
 * The state is the six doubles (Cg[0..5]) SSJ keeps; `set_seed` takes the same
 * six-long vector `setSeed(long[])` takes and applies the same validation, so a
 * call site can be transcribed from the Java verbatim.
 */
class Mrg32k3a {
public:
    Mrg32k3a() { reset_to_default(); }

    /** SSJ's setSeed(long[6]). Values are taken modulo the two moduli. */
    void set_seed(const long long s[6]) {
        validate(s);
        for (int i = 0; i < 6; ++i) cg_[i] = static_cast<double>(s[i]);
    }

    /** Convenience for the engine's `{seed+off, ..., seed+off+5}` idiom. */
    void set_seed_offset(long long seed, long long offset) {
        long long s[6];
        for (int i = 0; i < 6; ++i) s[i] = seed + offset + i;
        set_seed(s);
    }

    /**
     * SSJ's nextValue(): the combined generator, returning a double in (0, 1).
     * Never returns 0; may return values arbitrarily close to 1.
     */
    double next_double() {
        // Component 1
        double p1 = kA12 * cg_[1] - kA13n * cg_[0];
        double k = std::floor(p1 / kM1);
        p1 -= k * kM1;
        if (p1 < 0.0) p1 += kM1;
        cg_[0] = cg_[1];
        cg_[1] = cg_[2];
        cg_[2] = p1;

        // Component 2
        double p2 = kA21 * cg_[5] - kA23n * cg_[3];
        k = std::floor(p2 / kM2);
        p2 -= k * kM2;
        if (p2 < 0.0) p2 += kM2;
        cg_[3] = cg_[4];
        cg_[4] = cg_[5];
        cg_[5] = p2;

        // Combination
        return ((p1 > p2) ? (p1 - p2) * kNorm : (p1 - p2 + kM1) * kNorm);
    }

    /** SSJ's nextInt(i, j): a uniform integer on [i, j]. */
    int next_int(int i, int j) {
        if (i > j) throw InputError("Mrg32k3a::next_int: empty range");
        return i + static_cast<int>(next_double() * (static_cast<double>(j - i) + 1.0));
    }

    /** The six state doubles, for tests that pin the state and not only the draws. */
    double state(int i) const { return cg_[i]; }

private:
    void reset_to_default() {
        // SSJ's default initial seed (12345 in all six slots).
        for (int i = 0; i < 6; ++i) cg_[i] = 12345.0;
    }

    static void validate(const long long s[6]) {
        for (int i = 0; i < 3; ++i) {
            if (s[i] < 0 || static_cast<double>(s[i]) >= kM1)
                throw InputError("Mrg32k3a::set_seed: the first three seeds must be in [0, m1)");
        }
        for (int i = 3; i < 6; ++i) {
            if (s[i] < 0 || static_cast<double>(s[i]) >= kM2)
                throw InputError("Mrg32k3a::set_seed: the last three seeds must be in [0, m2)");
        }
        if (s[0] == 0 && s[1] == 0 && s[2] == 0)
            throw InputError("Mrg32k3a::set_seed: the first three seeds must not all be zero");
        if (s[3] == 0 && s[4] == 0 && s[5] == 0)
            throw InputError("Mrg32k3a::set_seed: the last three seeds must not all be zero");
    }

    // L'Ecuyer's constants, exactly as SSJ spells them.
    static constexpr double kM1 = 4294967087.0;
    static constexpr double kM2 = 4294944443.0;
    static constexpr double kA12 = 1403580.0;
    static constexpr double kA13n = 810728.0;
    static constexpr double kA21 = 527612.0;
    static constexpr double kA23n = 1370589.0;
    /** 1/(m1+1), SSJ's norm, so the result lies strictly inside (0,1). */
    static constexpr double kNorm = 2.328306549295727688e-10;

    double cg_[6];
};

// ---------------------------------------------------------------------------
// java.util.Random
// ---------------------------------------------------------------------------

/**
 * `java.util.Random`, the 48-bit LCG of the Java Language Specification.
 *
 * The engine uses 23 of these alongside the MRG streams, for the choices that
 * are draws rather than variates: which batch size a BMAP service releases,
 * which branch a probabilistic route takes, and so on. Transcribed rather than
 * approximated, for the same reason as above.
 */
class JavaRandom {
public:
    explicit JavaRandom(long long seed = 0) { set_seed(seed); }

    /** Java's setSeed: scramble by 0x5DEECE66D and mask to 48 bits. */
    void set_seed(long long seed) {
        seed_ = (static_cast<uint64_t>(seed) ^ 0x5DEECE66DULL) & ((1ULL << 48) - 1);
        have_next_gaussian_ = false;
    }

    /** Java's next(bits). */
    int32_t next(int bits) {
        seed_ = (seed_ * 0x5DEECE66DULL + 0xBULL) & ((1ULL << 48) - 1);
        return static_cast<int32_t>(static_cast<int64_t>(seed_) >> (48 - bits));
    }

    /** Java's nextInt(). */
    int32_t next_int() { return next(32); }

    /** Java's nextInt(bound), including the rejection loop it documents. */
    int32_t next_int(int32_t bound) {
        if (bound <= 0) throw InputError("JavaRandom::next_int: bound must be positive");
        if ((bound & -bound) == bound) {  // power of two
            return static_cast<int32_t>((static_cast<int64_t>(bound) * next(31)) >> 31);
        }
        int32_t bits, val;
        do {
            bits = next(31);
            val = bits % bound;
        } while (bits - val + (bound - 1) < 0);
        return val;
    }

    /** Java's nextDouble(): a 26-bit and a 27-bit draw. */
    double next_double() {
        return static_cast<double>((static_cast<int64_t>(next(26)) << 27) + next(27)) /
               static_cast<double>(1LL << 53);
    }

    /** Java's nextLong(). */
    int64_t next_long() {
        return (static_cast<int64_t>(next(32)) << 32) + next(32);
    }

    /** Java's nextBoolean(). */
    bool next_boolean() { return next(1) != 0; }

    /** Java's nextGaussian(), the polar method with its cached second value. */
    double next_gaussian() {
        if (have_next_gaussian_) {
            have_next_gaussian_ = false;
            return next_gaussian_;
        }
        double v1, v2, s;
        do {
            v1 = 2.0 * next_double() - 1.0;
            v2 = 2.0 * next_double() - 1.0;
            s = v1 * v1 + v2 * v2;
        } while (s >= 1.0 || s == 0.0);
        const double multiplier = std::sqrt(-2.0 * std::log(s) / s);
        next_gaussian_ = v2 * multiplier;
        have_next_gaussian_ = true;
        return v1 * multiplier;
    }

private:
    uint64_t seed_ = 0;
    bool have_next_gaussian_ = false;
    double next_gaussian_ = 0.0;
};

}  // namespace rng
}  // namespace line

#endif  // LINE_UTIL_RNG_SSJ_H
