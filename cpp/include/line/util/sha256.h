/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_SHA256_H
#define LINE_UTIL_SHA256_H

/**
 * @file
 * @ingroup line_util
 * SHA-256 (FIPS 180-4), for verifying a file the port did not build.
 *
 * The port needs this for the one thing it downloads, `JMT.jar`: the transfer
 * is delegated to curl or wget, which authenticate the SERVER, and a digest
 * over the received bytes is what authenticates the ARTEFACT. The two answer
 * different questions and the second is the one a jar about to be handed to a
 * JVM raises.
 *
 * Rather than link OpenSSL for one hash -- a link-time obligation on every
 * binary in the tree, for a dependency nothing else here wants -- this computes
 * it directly, as `util/websocket.h` already does for the SHA-1 of the
 * handshake. It is NOT a general cryptographic toolkit and must not become one:
 * the file digest below is its whole purpose.
 *
 * `sha256_file_hex` streams the file in blocks and never holds it in memory,
 * since the jar it was written for is 31 MiB.
 */

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>

namespace line {
namespace util {

namespace detail {

/** The 64 round constants: the cube roots of the first 64 primes, FIPS 180-4. */
static const std::uint32_t SHA256_K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u,
    0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu,
    0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu,
    0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau, 0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
    0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu,
    0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
    0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u, 0x19a4c116u,
    0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u,
    0xc67178f2u};

inline std::uint32_t sha256_ror(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

}  // namespace detail

/**
 * Incremental SHA-256. Feed it with update(), read the digest once with hex().
 *
 * The object is SPENT by hex(): the padding is appended to the running state
 * rather than to a copy, so a second call would digest the padding again. The
 * two free functions below are the intended interface; this class is public
 * only because a caller hashing something that is not a file or a string needs
 * somewhere to put the loop.
 */
class Sha256 {
public:
    Sha256() { reset(); }

    void reset() {
        h_[0] = 0x6a09e667u;
        h_[1] = 0xbb67ae85u;
        h_[2] = 0x3c6ef372u;
        h_[3] = 0xa54ff53au;
        h_[4] = 0x510e527fu;
        h_[5] = 0x9b05688cu;
        h_[6] = 0x1f83d9abu;
        h_[7] = 0x5be0cd19u;
        len_ = 0;
        total_ = 0;
    }

    void update(const void* data, std::size_t n) {
        const unsigned char* p = static_cast<const unsigned char*>(data);
        total_ += static_cast<std::uint64_t>(n);
        while (n > 0) {
            const std::size_t room = 64 - len_;
            const std::size_t take = n < room ? n : room;
            for (std::size_t i = 0; i < take; ++i) buf_[len_ + i] = p[i];
            len_ += take;
            p += take;
            n -= take;
            if (len_ == 64) {
                compress(buf_);
                len_ = 0;
            }
        }
    }

    void update(const std::string& s) { update(s.data(), s.size()); }

    /** Finalizes and returns the digest as 64 lowercase hex characters. */
    std::string hex() {
        const std::uint64_t bitlen = total_ * 8ull;
        unsigned char pad = 0x80;
        update(&pad, 1);
        pad = 0x00;
        while (len_ != 56) update(&pad, 1);
        unsigned char tail[8];
        for (int i = 0; i < 8; ++i)
            tail[i] = static_cast<unsigned char>((bitlen >> ((7 - i) * 8)) & 0xFFu);
        update(tail, 8);

        static const char* const digits = "0123456789abcdef";
        std::string out;
        out.reserve(64);
        for (int i = 0; i < 8; ++i)
            for (int b = 3; b >= 0; --b) {
                const unsigned char byte = static_cast<unsigned char>((h_[i] >> (b * 8)) & 0xFFu);
                out.push_back(digits[byte >> 4]);
                out.push_back(digits[byte & 0x0F]);
            }
        return out;
    }

private:
    void compress(const unsigned char block[64]) {
        std::uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = (std::uint32_t(block[i * 4]) << 24) | (std::uint32_t(block[i * 4 + 1]) << 16) |
                   (std::uint32_t(block[i * 4 + 2]) << 8) | std::uint32_t(block[i * 4 + 3]);
        for (int i = 16; i < 64; ++i) {
            const std::uint32_t s0 = detail::sha256_ror(w[i - 15], 7) ^
                                     detail::sha256_ror(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const std::uint32_t s1 = detail::sha256_ror(w[i - 2], 17) ^
                                     detail::sha256_ror(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        std::uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3];
        std::uint32_t e = h_[4], f = h_[5], g = h_[6], hh = h_[7];
        for (int i = 0; i < 64; ++i) {
            const std::uint32_t S1 =
                detail::sha256_ror(e, 6) ^ detail::sha256_ror(e, 11) ^ detail::sha256_ror(e, 25);
            const std::uint32_t ch = (e & f) ^ ((~e) & g);
            const std::uint32_t t1 = hh + S1 + ch + detail::SHA256_K[i] + w[i];
            const std::uint32_t S0 =
                detail::sha256_ror(a, 2) ^ detail::sha256_ror(a, 13) ^ detail::sha256_ror(a, 22);
            const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t t2 = S0 + maj;
            hh = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        h_[0] += a;
        h_[1] += b;
        h_[2] += c;
        h_[3] += d;
        h_[4] += e;
        h_[5] += f;
        h_[6] += g;
        h_[7] += hh;
    }

    std::uint32_t h_[8];
    unsigned char buf_[64];
    std::size_t len_;
    std::uint64_t total_;
};

/** The digest of a byte string, as 64 lowercase hex characters. */
inline std::string sha256_hex(const std::string& data) {
    Sha256 s;
    s.update(data);
    return s.hex();
}

/**
 * The digest of a file's contents, streamed.
 *
 * @return the 64 hex characters, or an EMPTY string when the file cannot be
 *         read -- which the caller must not confuse with a mismatch: a digest
 *         that could not be computed has verified nothing.
 */
inline std::string sha256_file_hex(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (f == nullptr) return std::string();
    Sha256 s;
    unsigned char buf[65536];
    while (true) {
        const std::size_t n = std::fread(buf, 1, sizeof(buf), f);
        if (n > 0) s.update(buf, n);
        if (n < sizeof(buf)) break;
    }
    const bool bad = std::ferror(f) != 0;
    std::fclose(f);
    return bad ? std::string() : s.hex();
}

}  // namespace util
}  // namespace line

#endif  // LINE_UTIL_SHA256_H
