/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_WEBSOCKET_H
#define LINE_UTIL_WEBSOCKET_H

/**
 * Minimal RFC 6455 WebSocket server, enough to serve `LineWebSocketServer`'s
 * protocol.
 *
 * The JAR runs `java -jar jline.jar -p <port>` as a solve server: a client
 * opens a WebSocket, sends ONE text message whose first line is the CSV
 * argument list and whose remainder is the model document, and receives the
 * CLI's output as one text message before the connection closes. This is the
 * same protocol, spoken directly over a POSIX socket.
 *
 * NO DEPENDENCY, for `http.h`'s reason: java-websocket is a link-time
 * obligation in the JAR and Boost.Beast would pull Asio and its thread
 * requirements into every binary in this tree, to speak a framing that fits in
 * two hundred lines. What is implemented is exactly what the protocol uses --
 * the opening handshake, masked client text frames of any length,
 * continuation frames, ping/pong and close.
 *
 * NO TLS, as `http.h` has none: the server binds a port for a local client, and
 * a caller who needs an authenticated channel must front it with a proxy rather
 * than believe this one provides it.
 *
 * ONE CONNECTION AT A TIME, deliberately. A solve is CPU-bound and the JAR's
 * server serves one message per connection and closes; accepting concurrently
 * would let two solves contend for the same cores and report timings neither
 * would produce alone.
 */

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

#include "line/util/error.h"

namespace line {
namespace ws {

/** Transport-level failure: the port is taken, the peer vanished, a bad frame. */
class WsError : public Error {
  public:
    explicit WsError(const std::string& what) : Error(what) {}
};

namespace detail {

/** SHA-1 of a byte string, RFC 3174; the handshake's only cryptographic need. */
inline std::string sha1(const std::string& msg) {
    std::uint32_t h[5] = {0x67452301u, 0xEFCDAB89u, 0x98BADCFEu, 0x10325476u, 0xC3D2E1F0u};
    std::string data = msg;
    const std::uint64_t bitlen = static_cast<std::uint64_t>(data.size()) * 8ull;
    data.push_back(static_cast<char>(0x80));
    while (data.size() % 64 != 56) data.push_back('\0');
    for (int i = 7; i >= 0; --i)
        data.push_back(static_cast<char>((bitlen >> (i * 8)) & 0xFF));

    for (std::size_t off = 0; off < data.size(); off += 64) {
        std::uint32_t w[80];
        for (int i = 0; i < 16; ++i) {
            const unsigned char* p =
                reinterpret_cast<const unsigned char*>(data.data() + off + i * 4);
            w[i] = (std::uint32_t(p[0]) << 24) | (std::uint32_t(p[1]) << 16) |
                   (std::uint32_t(p[2]) << 8) | std::uint32_t(p[3]);
        }
        for (int i = 16; i < 80; ++i) {
            const std::uint32_t v = w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16];
            w[i] = (v << 1) | (v >> 31);
        }
        std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4];
        for (int i = 0; i < 80; ++i) {
            std::uint32_t f = 0, k = 0;
            if (i < 20) {
                f = (b & c) | ((~b) & d);
                k = 0x5A827999u;
            } else if (i < 40) {
                f = b ^ c ^ d;
                k = 0x6ED9EBA1u;
            } else if (i < 60) {
                f = (b & c) | (b & d) | (c & d);
                k = 0x8F1BBCDCu;
            } else {
                f = b ^ c ^ d;
                k = 0xCA62C1D6u;
            }
            const std::uint32_t tmp = ((a << 5) | (a >> 27)) + f + e + k + w[i];
            e = d;
            d = c;
            c = (b << 30) | (b >> 2);
            b = a;
            a = tmp;
        }
        h[0] += a;
        h[1] += b;
        h[2] += c;
        h[3] += d;
        h[4] += e;
    }
    std::string out(20, '\0');
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 4; ++j)
            out[i * 4 + j] = static_cast<char>((h[i] >> ((3 - j) * 8)) & 0xFF);
    return out;
}

inline std::string base64(const std::string& in) {
    static const char* tbl = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    std::size_t i = 0;
    while (i + 2 < in.size()) {
        const unsigned v = (static_cast<unsigned char>(in[i]) << 16) |
                           (static_cast<unsigned char>(in[i + 1]) << 8) |
                           static_cast<unsigned char>(in[i + 2]);
        out.push_back(tbl[(v >> 18) & 63]);
        out.push_back(tbl[(v >> 12) & 63]);
        out.push_back(tbl[(v >> 6) & 63]);
        out.push_back(tbl[v & 63]);
        i += 3;
    }
    if (i + 1 == in.size()) {
        const unsigned v = static_cast<unsigned char>(in[i]) << 16;
        out.push_back(tbl[(v >> 18) & 63]);
        out.push_back(tbl[(v >> 12) & 63]);
        out += "==";
    } else if (i + 2 == in.size()) {
        const unsigned v = (static_cast<unsigned char>(in[i]) << 16) |
                           (static_cast<unsigned char>(in[i + 1]) << 8);
        out.push_back(tbl[(v >> 18) & 63]);
        out.push_back(tbl[(v >> 12) & 63]);
        out.push_back(tbl[(v >> 6) & 63]);
        out.push_back('=');
    }
    return out;
}

inline bool send_all(int fd, const char* p, std::size_t n) {
    while (n) {
        const ssize_t k = ::send(fd, p, n, MSG_NOSIGNAL);
        if (k <= 0) return false;
        p += k;
        n -= static_cast<std::size_t>(k);
    }
    return true;
}

inline bool recv_exact(int fd, char* p, std::size_t n) {
    while (n) {
        const ssize_t k = ::recv(fd, p, n, 0);
        if (k <= 0) return false;
        p += k;
        n -= static_cast<std::size_t>(k);
    }
    return true;
}

/** Case-insensitive header lookup over a raw request head. */
inline std::string header(const std::string& head, const std::string& key) {
    std::string lower = head, lkey = key;
    for (std::size_t i = 0; i < lower.size(); ++i)
        lower[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(lower[i])));
    for (std::size_t i = 0; i < lkey.size(); ++i)
        lkey[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(lkey[i])));
    std::size_t at = lower.find("\n" + lkey + ":");
    if (at == std::string::npos) return std::string();
    at += lkey.size() + 2;
    const std::size_t end = head.find('\n', at);
    std::string v = head.substr(at, end == std::string::npos ? std::string::npos : end - at);
    while (!v.empty() && (v[0] == ' ' || v[0] == '\t')) v.erase(v.begin());
    while (!v.empty() && (v.back() == '\r' || v.back() == ' ')) v.pop_back();
    return v;
}

}  // namespace detail

/** A listening socket; one connection is served at a time. */
class Server {
  public:
    explicit Server(int port) {
        fd_ = ::socket(AF_INET, SOCK_STREAM, 0);
        if (fd_ < 0) throw WsError("cannot create a listening socket");
        int on = 1;
        ::setsockopt(fd_, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on));
        sockaddr_in a;
        std::memset(&a, 0, sizeof(a));
        a.sin_family = AF_INET;
        a.sin_addr.s_addr = htonl(INADDR_ANY);
        a.sin_port = htons(static_cast<uint16_t>(port));
        if (::bind(fd_, reinterpret_cast<sockaddr*>(&a), sizeof(a)) != 0) {
            ::close(fd_);
            throw WsError("cannot bind port " + std::to_string(port) +
                          "; another process is probably listening on it");
        }
        if (::listen(fd_, 4) != 0) {
            ::close(fd_);
            throw WsError("cannot listen on port " + std::to_string(port));
        }
    }
    ~Server() {
        if (fd_ >= 0) ::close(fd_);
    }
    Server(const Server&) = delete;
    Server& operator=(const Server&) = delete;

    /**
     * Accept one connection, complete the handshake, read ONE text message and
     * hand it to `serve`; send what `serve` returns and close.
     *
     * @return false when the peer failed before a message arrived, which is a
     *         dropped client and not a reason to stop the server.
     */
    template <class Fn>
    bool serve_one(Fn serve) {
        const int c = ::accept(fd_, nullptr, nullptr);
        if (c < 0) return false;
        const bool ok = handshake(c) && exchange(c, serve);
        ::close(c);
        return ok;
    }

  private:
    int fd_ = -1;

    static bool handshake(int c) {
        std::string head;
        char buf[1024];
        while (head.find("\r\n\r\n") == std::string::npos) {
            const ssize_t k = ::recv(c, buf, sizeof(buf), 0);
            if (k <= 0) return false;
            head.append(buf, static_cast<std::size_t>(k));
            if (head.size() > 65536) return false;  // not a handshake
        }
        const std::string key = detail::header(head, "Sec-WebSocket-Key");
        if (key.empty()) return false;
        // The RFC's fixed GUID; the accept token is base64(sha1(key + GUID)).
        const std::string accept = detail::base64(
            detail::sha1(key + "258EAFA5-E914-47DA-95CA-C5AB0DC85B11"));
        const std::string resp = "HTTP/1.1 101 Switching Protocols\r\n"
                                 "Upgrade: websocket\r\n"
                                 "Connection: Upgrade\r\n"
                                 "Sec-WebSocket-Accept: " + accept + "\r\n\r\n";
        return detail::send_all(c, resp.data(), resp.size());
    }

    /** Read one complete text message, possibly fragmented. */
    static bool read_message(int c, std::string& out) {
        out.clear();
        for (;;) {
            unsigned char h[2];
            if (!detail::recv_exact(c, reinterpret_cast<char*>(h), 2)) return false;
            const bool fin = (h[0] & 0x80) != 0;
            const int opcode = h[0] & 0x0F;
            const bool masked = (h[1] & 0x80) != 0;
            std::uint64_t len = h[1] & 0x7F;
            if (len == 126) {
                unsigned char e[2];
                if (!detail::recv_exact(c, reinterpret_cast<char*>(e), 2)) return false;
                len = (std::uint64_t(e[0]) << 8) | e[1];
            } else if (len == 127) {
                unsigned char e[8];
                if (!detail::recv_exact(c, reinterpret_cast<char*>(e), 8)) return false;
                len = 0;
                for (int i = 0; i < 8; ++i) len = (len << 8) | e[i];
            }
            // A client frame MUST be masked (RFC 6455 5.1); an unmasked one is
            // a protocol error, not a frame to decode as if the mask were zero.
            unsigned char mask[4] = {0, 0, 0, 0};
            if (masked && !detail::recv_exact(c, reinterpret_cast<char*>(mask), 4)) return false;
            if (!masked) return false;
            if (len > (1ull << 30)) return false;  // a model document, not a stream
            std::string payload(static_cast<std::size_t>(len), '\0');
            if (len && !detail::recv_exact(c, &payload[0], payload.size())) return false;
            for (std::size_t i = 0; i < payload.size(); ++i)
                payload[i] = static_cast<char>(payload[i] ^ mask[i % 4]);

            if (opcode == 0x8) return false;  // close
            if (opcode == 0x9) {              // ping -> pong, then keep reading
                std::string pong;
                pong.push_back(static_cast<char>(0x8A));
                pong.push_back(static_cast<char>(payload.size() & 0x7F));
                pong += payload;
                if (!detail::send_all(c, pong.data(), pong.size())) return false;
                continue;
            }
            if (opcode == 0xA) continue;  // pong
            out += payload;
            if (fin) return true;
        }
    }

    static bool send_text(int c, const std::string& msg) {
        std::string f;
        f.push_back(static_cast<char>(0x81));  // FIN + text
        if (msg.size() < 126) {
            f.push_back(static_cast<char>(msg.size()));
        } else if (msg.size() <= 0xFFFF) {
            f.push_back(static_cast<char>(126));
            f.push_back(static_cast<char>((msg.size() >> 8) & 0xFF));
            f.push_back(static_cast<char>(msg.size() & 0xFF));
        } else {
            f.push_back(static_cast<char>(127));
            for (int i = 7; i >= 0; --i)
                f.push_back(static_cast<char>((static_cast<std::uint64_t>(msg.size()) >> (i * 8)) &
                                              0xFF));
        }
        f += msg;
        if (!detail::send_all(c, f.data(), f.size())) return false;
        // A clean close, so the client sees the message as complete rather than
        // as a connection that dropped mid-answer.
        const char close_frame[2] = {static_cast<char>(0x88), 0};
        detail::send_all(c, close_frame, 2);
        return true;
    }

    template <class Fn>
    static bool exchange(int c, Fn serve) {
        std::string msg;
        if (!read_message(c, msg)) return false;
        return send_text(c, serve(msg));
    }
};

}  // namespace ws
}  // namespace line

#endif  // LINE_UTIL_WEBSOCKET_H
