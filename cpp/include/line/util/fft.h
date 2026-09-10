/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_FFT_H
#define LINE_UTIL_FFT_H

/**
 * Discrete Fourier transform of arbitrary length, in double complex.
 *
 * WHY IT IS HERE. `MG1_CR`, the Bini-Meini point-wise cyclic reduction behind
 * every M/G/1-type G matrix, evaluates four matrix polynomials at the (nj+1)-th
 * roots of unity, combines them through a point-wise inverse, and interpolates
 * the result back. The transform length is nj+1 where nj+1 doubles from the
 * degree of the block sequence, so it is NOT a power of two in general: a
 * BMAP with three batch sizes gives a degree-4 sequence and a first transform
 * of length 4, but a degree-5 one gives length 6. Both must work, and both must
 * agree with MATLAB's `fft`/`ifft` to roundoff, because the reference's
 * truncation test compares the interpolated tail against an absolute epsilon.
 *
 * The implementation is radix-2 Cooley-Tukey when the length is a power of two
 * and Bluestein's chirp-z otherwise, which reduces the arbitrary length to a
 * power-of-two convolution. Bluestein is exact in the same sense the radix-2
 * transform is -- it is a rearrangement, not an approximation -- so the two
 * paths differ only in rounding, and neither introduces the O(n) error a naive
 * quadratic DFT accumulates at n = 2048 (the reference's MaxNumRoot).
 *
 * CONVENTION, matching MATLAB. `dft(a, false)` returns
 * sum_k a_k exp(-2 pi i j k / n); `dft(a, true)` returns
 * (1/n) sum_k a_k exp(+2 pi i j k / n). The inverse therefore carries the 1/n,
 * exactly as `ifft` does, so a forward followed by an inverse is the identity.
 *
 * DOUBLE ONLY. The twiddle factors are cos/sin, which the exact and
 * multiprecision number types do not provide; callers that need the transform
 * at another arithmetic must say so rather than silently losing precision here.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "line/util/error.h"

namespace line {

namespace fft_detail {

/** Iterative radix-2 Cooley-Tukey, in place. Length must be a power of two. */
inline void radix2(std::vector<std::complex<double>>& a, bool conjugate) {
    const std::size_t n = a.size();
    if (n <= 1) return;
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    const double pi = 3.14159265358979323846;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2.0 * pi / static_cast<double>(len) * (conjugate ? 1.0 : -1.0);
        const std::complex<double> wlen(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (std::size_t k = 0; k < len / 2; ++k) {
                const std::complex<double> u = a[i + k];
                const std::complex<double> v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

inline bool is_power_of_two(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

/** Bluestein's chirp-z transform for a length that is not a power of two. */
inline void bluestein(std::vector<std::complex<double>>& a, bool conjugate) {
    const std::size_t n = a.size();
    const double pi = 3.14159265358979323846;
    const double sign = conjugate ? 1.0 : -1.0;

    // The chirp exp(sign i pi k^2 / n); k^2 is reduced mod 2n so the angle
    // stays small and the cos/sin keep their relative accuracy at large k.
    std::vector<std::complex<double>> chirp(n);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t kk = (k * k) % (2 * n);
        const double ang = sign * pi * static_cast<double>(kk) / static_cast<double>(n);
        chirp[k] = std::complex<double>(std::cos(ang), std::sin(ang));
    }

    std::size_t m = 1;
    while (m < 2 * n - 1) m <<= 1;
    std::vector<std::complex<double>> x(m, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> y(m, std::complex<double>(0.0, 0.0));
    for (std::size_t k = 0; k < n; ++k) x[k] = a[k] * chirp[k];
    y[0] = std::conj(chirp[0]);
    for (std::size_t k = 1; k < n; ++k) {
        y[k] = std::conj(chirp[k]);
        y[m - k] = std::conj(chirp[k]);
    }

    radix2(x, false);
    radix2(y, false);
    for (std::size_t k = 0; k < m; ++k) x[k] *= y[k];
    radix2(x, true);
    const double scale = 1.0 / static_cast<double>(m);
    for (std::size_t k = 0; k < n; ++k) a[k] = x[k] * scale * chirp[k];
}

}  // namespace fft_detail

/**
 * In-place DFT of `a`. `inverse` selects the +i sign and the 1/n scaling, so
 * the pair matches MATLAB's `fft` and `ifft`.
 */
inline void dft(std::vector<std::complex<double>>& a, bool inverse) {
    const std::size_t n = a.size();
    if (n <= 1) return;
    if (fft_detail::is_power_of_two(n)) {
        fft_detail::radix2(a, inverse);
    } else {
        fft_detail::bluestein(a, inverse);
    }
    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (std::size_t k = 0; k < n; ++k) a[k] *= scale;
    }
}

}  // namespace line

#endif  // LINE_UTIL_FFT_H
