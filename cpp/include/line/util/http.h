/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_HTTP_H
#define LINE_UTIL_HTTP_H

/**
 * Minimal HTTP/1.1 client, enough to talk to a line-*-rest service.
 *
 * The port needs this for api/sym, whose only backend is the line-sage-rest
 * service: the JAR reaches it through HttpURLConnection, MATLAB through
 * webwrite and Python through requests, and there is no equivalent in the C++
 * standard library. Rather than take a dependency on libcurl or Boost.Beast --
 * the first is a link-time obligation for every binary in the tree, the second
 * pulls in Asio and its thread requirements -- this speaks the protocol
 * directly over a POSIX socket. What it supports is exactly what the service
 * uses: GET, POST of a JSON body, identity and chunked transfer encodings,
 * connect and read timeouts.
 *
 * NO TLS. https:// is REFUSED BY NAME rather than silently downgraded to http.
 * That matters because the symbolic backend is normally a container on
 * localhost, where plaintext is correct, and a caller pointing LINE_SAGE_URL at
 * a remote https endpoint must be told that this client cannot authenticate it
 * rather than have the request fail obscurely or, worse, travel in the clear.
 *
 * Redirects are NOT followed and cookies are not kept: the REST protocol has
 * neither, and following a redirect silently would let a 302 turn a POST into a
 * GET and lose the request body.
 */

#include <cctype>
#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>

#include <fcntl.h>
#include <netdb.h>
#include <poll.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

#include "line/util/error.h"

namespace line {
namespace http {

/** Transport-level failure: unresolvable host, refused connection, timeout. */
class HttpError : public Error {
public:
    explicit HttpError(const std::string& what) : Error(what) {}
};

/** An HTTP response, with the body already de-chunked. */
struct Response {
    int status = 0;    ///< HTTP status code
    std::string body;  ///< Response body, decoded
};

/** The pieces of an http:// URL this client needs. */
struct Url {
    std::string host;
    int port = 80;
    std::string path = "/";
};

/**
 * Splits an http:// URL. Supports host, host:port and [v6addr]:port forms.
 *
 * @param url the URL
 * @return its host, port and path
 */
inline Url parse_url(const std::string& url) {
    if (url.compare(0, 8, "https://") == 0)
        throw UnsupportedError(
            "line::http: https is not supported by this client, which has no TLS; point the "
            "service URL at an http:// endpoint (a local container is the usual case) or "
            "terminate TLS in front of it");
    if (url.compare(0, 7, "http://") != 0)
        throw InputError("line::http: the URL must start with http://, got '" + url + "'");

    const std::string rest = url.substr(7);
    const std::size_t slash = rest.find('/');
    const std::string authority = slash == std::string::npos ? rest : rest.substr(0, slash);
    Url u;
    u.path = slash == std::string::npos ? "/" : rest.substr(slash);
    if (u.path.empty()) u.path = "/";

    if (!authority.empty() && authority[0] == '[') {
        const std::size_t close = authority.find(']');
        if (close == std::string::npos)
            throw InputError("line::http: unterminated IPv6 literal in '" + url + "'");
        u.host = authority.substr(1, close - 1);
        if (close + 1 < authority.size() && authority[close + 1] == ':')
            u.port = std::atoi(authority.c_str() + close + 2);
    } else {
        const std::size_t colon = authority.rfind(':');
        if (colon == std::string::npos) {
            u.host = authority;
        } else {
            u.host = authority.substr(0, colon);
            u.port = std::atoi(authority.c_str() + colon + 1);
        }
    }
    if (u.host.empty()) throw InputError("line::http: no host in '" + url + "'");
    if (u.port <= 0 || u.port > 65535)
        throw InputError("line::http: port out of range in '" + url + "'");
    return u;
}

namespace detail {

/** RAII holder so every early return closes the socket. */
class Socket {
public:
    explicit Socket(int fd = -1) : fd_(fd) {}
    ~Socket() {
        if (fd_ >= 0) ::close(fd_);
    }
    Socket(const Socket&) = delete;
    Socket& operator=(const Socket&) = delete;
    int fd() const { return fd_; }

private:
    int fd_;
};

/**
 * Connects with a bounded wait. The connect timeout is enforced by hand through
 * a nonblocking connect plus poll(): SO_SNDTIMEO does NOT bound connect() on
 * Linux, so relying on it would let an unreachable host hang for the kernel's
 * own SYN retry budget, over two minutes.
 */
inline int connect_socket(const std::string& host, int port, int connectMillis) {
    struct addrinfo hints;
    std::memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;

    struct addrinfo* res = nullptr;
    const std::string portStr = std::to_string(port);
    if (::getaddrinfo(host.c_str(), portStr.c_str(), &hints, &res) != 0 || res == nullptr)
        throw HttpError("line::http: cannot resolve " + host + ":" + portStr);

    int connected = -1;
    std::string lastError = "connection refused";
    for (struct addrinfo* ai = res; ai != nullptr && connected < 0; ai = ai->ai_next) {
        const int fd = ::socket(ai->ai_family, ai->ai_socktype, ai->ai_protocol);
        if (fd < 0) continue;
        const int flags = ::fcntl(fd, F_GETFL, 0);
        ::fcntl(fd, F_SETFL, flags | O_NONBLOCK);

        int rc = ::connect(fd, ai->ai_addr, ai->ai_addrlen);
        if (rc != 0 && errno == EINPROGRESS) {
            struct pollfd pfd;
            pfd.fd = fd;
            pfd.events = POLLOUT;
            pfd.revents = 0;
            const int pr = ::poll(&pfd, 1, connectMillis > 0 ? connectMillis : -1);
            if (pr > 0) {
                int soerr = 0;
                socklen_t len = sizeof(soerr);
                if (::getsockopt(fd, SOL_SOCKET, SO_ERROR, &soerr, &len) == 0 && soerr == 0)
                    rc = 0;
                else
                    lastError = std::strerror(soerr);
            } else if (pr == 0) {
                lastError = "connect timed out";
                rc = -1;
            } else {
                lastError = std::strerror(errno);
                rc = -1;
            }
        } else if (rc != 0) {
            lastError = std::strerror(errno);
        }

        if (rc == 0) {
            ::fcntl(fd, F_SETFL, flags);
            connected = fd;
        } else {
            ::close(fd);
        }
    }
    ::freeaddrinfo(res);
    if (connected < 0)
        throw HttpError("line::http: cannot connect to " + host + ":" + portStr + " (" +
                        lastError + ")");
    return connected;
}

/** Applies a receive and send timeout to an established socket. */
inline void set_io_timeout(int fd, int millis) {
    if (millis <= 0) return;
    struct timeval tv;
    tv.tv_sec = millis / 1000;
    tv.tv_usec = (millis % 1000) * 1000;
    ::setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    ::setsockopt(fd, SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
}

inline void send_all(int fd, const std::string& data) {
    std::size_t sent = 0;
    while (sent < data.size()) {
#ifdef MSG_NOSIGNAL
        const ssize_t k = ::send(fd, data.data() + sent, data.size() - sent, MSG_NOSIGNAL);
#else
        const ssize_t k = ::send(fd, data.data() + sent, data.size() - sent, 0);
#endif
        if (k > 0) {
            sent += static_cast<std::size_t>(k);
            continue;
        }
        if (k < 0 && (errno == EINTR)) continue;
        throw HttpError(std::string("line::http: write failed (") + std::strerror(errno) + ")");
    }
}

/** Appends one read to buf. Returns false at end of stream. */
inline bool read_more(int fd, std::string& buf) {
    char chunk[8192];
    while (true) {
        const ssize_t k = ::recv(fd, chunk, sizeof(chunk), 0);
        if (k > 0) {
            buf.append(chunk, static_cast<std::size_t>(k));
            return true;
        }
        if (k == 0) return false;
        if (errno == EINTR) continue;
        if (errno == EAGAIN || errno == EWOULDBLOCK)
            throw HttpError("line::http: read timed out");
        throw HttpError(std::string("line::http: read failed (") + std::strerror(errno) + ")");
    }
}

/** Case-insensitive header lookup over the raw header block. */
inline std::string header_value(const std::string& headers, const std::string& name) {
    std::string lowerHeaders(headers);
    for (std::size_t i = 0; i < lowerHeaders.size(); ++i)
        lowerHeaders[i] = static_cast<char>(std::tolower(lowerHeaders[i]));
    std::string key = "\r\n" + name + ":";
    for (std::size_t i = 0; i < key.size(); ++i)
        key[i] = static_cast<char>(std::tolower(key[i]));
    const std::size_t at = lowerHeaders.find(key);
    if (at == std::string::npos) return std::string();
    const std::size_t from = at + key.size();
    const std::size_t eol = headers.find("\r\n", from);
    std::string v = headers.substr(from, eol == std::string::npos ? std::string::npos : eol - from);
    const std::size_t b = v.find_first_not_of(" \t");
    const std::size_t e = v.find_last_not_of(" \t");
    return b == std::string::npos ? std::string() : v.substr(b, e - b + 1);
}

/** Reads the status line, the headers and the body, de-chunking if needed. */
inline Response read_response(int fd) {
    std::string buf;
    std::size_t headerEnd = std::string::npos;
    while ((headerEnd = buf.find("\r\n\r\n")) == std::string::npos) {
        if (!read_more(fd, buf)) break;
    }
    if (headerEnd == std::string::npos)
        throw HttpError("line::http: the server closed the connection before sending headers");

    // the leading CRLF lets header_value() match the first header too
    const std::string headers = "\r\n" + buf.substr(0, headerEnd + 2);
    std::string body = buf.substr(headerEnd + 4);

    Response r;
    const std::size_t sp = buf.find(' ');
    if (sp == std::string::npos) throw HttpError("line::http: malformed status line");
    r.status = std::atoi(buf.c_str() + sp + 1);

    const std::string encoding = header_value(headers, "transfer-encoding");
    const std::string length = header_value(headers, "content-length");
    if (encoding.find("chunked") != std::string::npos) {
        std::string decoded;
        std::size_t at = 0;
        while (true) {
            std::size_t eol;
            while ((eol = body.find("\r\n", at)) == std::string::npos) {
                if (!read_more(fd, body))
                    throw HttpError("line::http: truncated chunked body");
            }
            const std::size_t size =
                static_cast<std::size_t>(std::strtoul(body.c_str() + at, nullptr, 16));
            at = eol + 2;
            if (size == 0) break;
            while (body.size() < at + size + 2) {
                if (!read_more(fd, body))
                    throw HttpError("line::http: truncated chunked body");
            }
            decoded.append(body, at, size);
            at += size + 2;
        }
        r.body = decoded;
    } else if (!length.empty()) {
        const std::size_t want = static_cast<std::size_t>(std::strtoul(length.c_str(), nullptr, 10));
        while (body.size() < want) {
            if (!read_more(fd, body)) break;
        }
        r.body = body.substr(0, want < body.size() ? want : body.size());
    } else {
        while (read_more(fd, body)) {
        }
        r.body = body;
    }
    return r;
}

/** Sends one request on a fresh connection and reads the whole response. */
inline Response request(const std::string& method, const std::string& url, const std::string& body,
                        const std::string& contentType, int timeoutMillis) {
    const Url u = parse_url(url);
    const int connectMillis =
        timeoutMillis > 0 && timeoutMillis < 30000 ? timeoutMillis : 30000;
    Socket sock(connect_socket(u.host, u.port, connectMillis));
    set_io_timeout(sock.fd(), timeoutMillis);

    std::string req = method + " " + u.path + " HTTP/1.1\r\n";
    req += "Host: " + u.host + ":" + std::to_string(u.port) + "\r\n";
    req += "User-Agent: line-cpp\r\n";
    req += "Accept: application/json\r\n";
    req += "Connection: close\r\n";
    if (!body.empty()) {
        req += "Content-Type: " + contentType + "\r\n";
        req += "Content-Length: " + std::to_string(body.size()) + "\r\n";
    }
    req += "\r\n";
    req += body;

    send_all(sock.fd(), req);
    return read_response(sock.fd());
}

}  // namespace detail

/**
 * GET a URL.
 *
 * @param url           the absolute http:// URL
 * @param timeoutMillis read timeout in milliseconds, 0 for none
 * @return the response
 */
inline Response get(const std::string& url, int timeoutMillis) {
    return detail::request("GET", url, std::string(), std::string(), timeoutMillis);
}

/**
 * POST a JSON document.
 *
 * @param url           the absolute http:// URL
 * @param json          the request body
 * @param timeoutMillis read timeout in milliseconds, 0 for none
 * @return the response
 */
inline Response post_json(const std::string& url, const std::string& json, int timeoutMillis) {
    return detail::request("POST", url, json, "application/json; charset=utf-8", timeoutMillis);
}

}  // namespace http
}  // namespace line

#endif  // LINE_UTIL_HTTP_H
