/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AG_AG_WORKER_CLIENT_H
#define LINE_SOLVERS_AG_AG_WORKER_CLIENT_H

/**
 * @file ag_worker_client.h
 * @brief Coordinator-side connections to the ag-worker processes.
 *
 * Speaks the protocol of `jline.solvers.ag.AgWire`: newline-delimited JSON over
 * a plain TCP socket, one object per line. `assign` ships the STATIC half of
 * each owned agent -- its local rate matrix and the passive/active matrices of
 * the actions it takes part in -- once per solve, because none of it changes
 * across the fixed point; `sweep` then carries only the reversed rates, one
 * double per action, which is the entire coupling between agents.
 *
 * A worker that cannot be reached, or that fails mid-run, is marked dead and its
 * agents fall back to the coordinator. That is not a defensive workaround but a
 * property of the decomposition: an agent depends on the rest of the model only
 * through the reversed rates, so anyone holding them can solve it.
 */

#include <array>
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <unistd.h>

#include "json.hpp"

#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ag {

/** One agent's static description, in the wire form the worker expects. */
struct AgWirePayload {
    int k = 0;
    int n = 0;
    int mph = 1;
    int nlev = 1;
    std::vector<int> level;
    /** [row, col, value] triplets, 0-based. */
    std::vector<std::array<double, 3>> L;
    std::vector<int> passive_c;
    std::vector<std::vector<std::array<double, 3>>> passive_m;
    std::vector<int> active_c;
    std::vector<std::vector<std::array<double, 3>>> active_m;
};

/** Non-zero entries of M as [row, col, value] triplets, 0-based. */
inline std::vector<std::array<double, 3>> ag_triplets(const Matrix<double>& m) {
    std::vector<std::array<double, 3>> out;
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            const double v = m(i, j);
            if (v != 0.0) {
                out.push_back({static_cast<double>(i), static_cast<double>(j), v});
            }
        }
    }
    return out;
}

/** The connections to every configured worker, plus the agent partition. */
class AgWorkerPool {
  public:
    AgWorkerPool(const std::vector<std::string>& endpoints, double timeout_seconds)
        : endpoints_(endpoints), timeout_(timeout_seconds) {
        if (endpoints_.empty()) {
            throw InputError(
                "ag: the 'cluster' execution backend needs worker endpoints: set "
                "AgOptions::endpoints to \"host:port\" strings, each one an ag-worker "
                "started with 'java -cp jline.jar jline.cli.AgWorker -p <port>'");
        }
        fds_.assign(endpoints_.size(), -1);
        live_.assign(endpoints_.size(), true);
    }

    ~AgWorkerPool() { close_all(); }

    AgWorkerPool(const AgWorkerPool&) = delete;
    AgWorkerPool& operator=(const AgWorkerPool&) = delete;

    /**
     * Connect and ship the static half of each owned agent, once. @p payload maps
     * an agent index to its wire description.
     */
    template <class Payload>
    void ensure_assigned(std::size_t num_agents, Payload payload) {
        if (assigned_) return;
        partition(num_agents);
        for (std::size_t w = 0; w < endpoints_.size(); ++w) {
            if (owns_[w].empty()) continue;
            try {
                connect(w);
                nlohmann::json agents = nlohmann::json::array();
                for (std::size_t t = 0; t < owns_[w].size(); ++t) {
                    agents.push_back(encode(payload(owns_[w][t])));
                }
                nlohmann::json msg;
                msg["op"] = "assign";
                msg["agents"] = agents;
                const nlohmann::json reply = call(w, msg);
                if (reply.value("op", std::string()) != "assigned") {
                    throw std::runtime_error("worker did not acknowledge the assignment");
                }
            } catch (const std::exception& e) {
                std::cerr << "[LINE] Warning: AG worker " << endpoints_[w]
                          << " is unreachable (" << e.what()
                          << "); its agents run on the coordinator instead" << std::endl;
                live_[w] = false;
                close_one(w);
            }
        }
        assigned_ = true;
    }

    /**
     * One sweep. Fills @p pi for every agent a live worker answered for and
     * clears its @p pending flag; the rest are left to the caller.
     */
    void sweep(const std::vector<double>& x,
               std::vector<std::vector<double>>& pi,
               std::vector<bool>& pending) {
        for (std::size_t w = 0; w < endpoints_.size(); ++w) {
            if (!live_[w] || owns_[w].empty()) continue;
            try {
                nlohmann::json msg;
                msg["op"] = "sweep";
                msg["x"] = x;
                const nlohmann::json reply = call(w, msg);
                if (reply.value("op", std::string()) != "swept") {
                    throw std::runtime_error("unexpected reply '" +
                                             reply.value("op", std::string()) + "'");
                }
                for (const nlohmann::json& e : reply.at("agents")) {
                    const std::size_t k = e.at("k").get<std::size_t>();
                    if (k >= pi.size()) {
                        throw std::runtime_error("worker answered for an agent out of range");
                    }
                    pi[k] = e.at("pi").get<std::vector<double>>();
                    pending[k] = false;
                }
            } catch (const std::exception& e) {
                std::cerr << "[LINE] Warning: AG worker " << endpoints_[w]
                          << " failed mid-sweep (" << e.what()
                          << "); its agents are solved on the coordinator for the rest of "
                          << "the run" << std::endl;
                live_[w] = false;
                close_one(w);
            }
        }
    }

  private:
    /**
     * Round-robin in agent index order, computed before any connection is
     * attempted so a dead worker does not shift the others' agents. A pure
     * function of (agent count, worker count), so a rerun assigns the same
     * agents to the same workers.
     */
    void partition(std::size_t n) {
        owns_.assign(endpoints_.size(), std::vector<std::size_t>());
        for (std::size_t k = 0; k < n; ++k) owns_[k % endpoints_.size()].push_back(k);
    }

    static nlohmann::json encode(const AgWirePayload& a) {
        nlohmann::json j;
        j["k"] = a.k;
        j["n"] = a.n;
        j["mph"] = a.mph;
        j["nlev"] = a.nlev;
        j["level"] = a.level;
        j["L"] = a.L;
        nlohmann::json passive = nlohmann::json::array();
        for (std::size_t i = 0; i < a.passive_c.size(); ++i) {
            nlohmann::json e;
            e["c"] = a.passive_c[i];
            e["M"] = a.passive_m[i];
            passive.push_back(e);
        }
        nlohmann::json active = nlohmann::json::array();
        for (std::size_t i = 0; i < a.active_c.size(); ++i) {
            nlohmann::json e;
            e["c"] = a.active_c[i];
            e["M"] = a.active_m[i];
            active.push_back(e);
        }
        j["passive"] = passive;
        j["active"] = active;
        return j;
    }

    void connect(std::size_t w) {
        const std::string& ep = endpoints_[w];
        const std::size_t colon = ep.rfind(':');
        if (colon == std::string::npos || colon == 0) {
            throw std::runtime_error("malformed endpoint '" + ep + "', expected host:port");
        }
        const std::string host = ep.substr(0, colon);
        const std::string port = ep.substr(colon + 1);

        addrinfo hints{};
        hints.ai_family = AF_UNSPEC;
        hints.ai_socktype = SOCK_STREAM;
        addrinfo* res = nullptr;
        if (::getaddrinfo(host.c_str(), port.c_str(), &hints, &res) != 0 || res == nullptr) {
            throw std::runtime_error("cannot resolve " + ep);
        }
        int fd = -1;
        for (addrinfo* p = res; p != nullptr; p = p->ai_next) {
            fd = ::socket(p->ai_family, p->ai_socktype, p->ai_protocol);
            if (fd < 0) continue;
            timeval tv{};
            tv.tv_sec = static_cast<long>(timeout_);
            tv.tv_usec = static_cast<long>((timeout_ - tv.tv_sec) * 1e6);
            ::setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
            ::setsockopt(fd, SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
            if (::connect(fd, p->ai_addr, p->ai_addrlen) == 0) break;
            ::close(fd);
            fd = -1;
        }
        ::freeaddrinfo(res);
        if (fd < 0) throw std::runtime_error("cannot connect to " + ep);
        fds_[w] = fd;
    }

    nlohmann::json call(std::size_t w, const nlohmann::json& msg) {
        const std::string line = msg.dump() + "\n";
        std::size_t sent = 0;
        while (sent < line.size()) {
            const ssize_t n = ::send(fds_[w], line.data() + sent, line.size() - sent, 0);
            if (n <= 0) throw std::runtime_error("send failed: " + std::string(std::strerror(errno)));
            sent += static_cast<std::size_t>(n);
        }
        std::string reply;
        char ch = 0;
        while (true) {
            const ssize_t n = ::recv(fds_[w], &ch, 1, 0);
            if (n <= 0) throw std::runtime_error("worker closed the connection");
            if (ch == '\n') break;
            reply.push_back(ch);
        }
        nlohmann::json j = nlohmann::json::parse(reply);
        if (j.contains("error")) throw std::runtime_error(j.at("error").get<std::string>());
        return j;
    }

    void close_one(std::size_t w) {
        if (fds_[w] >= 0) {
            ::close(fds_[w]);
            fds_[w] = -1;
        }
    }

    void close_all() {
        for (std::size_t w = 0; w < fds_.size(); ++w) {
            if (fds_[w] < 0) continue;
            const std::string bye = "{\"op\":\"bye\"}\n";
            ::send(fds_[w], bye.data(), bye.size(), 0);
            close_one(w);
        }
    }

    std::vector<std::string> endpoints_;
    double timeout_;
    std::vector<int> fds_;
    std::vector<bool> live_;
    std::vector<std::vector<std::size_t>> owns_;
    bool assigned_ = false;
};

}  // namespace ag
}  // namespace line

#endif  // LINE_SOLVERS_AG_AG_WORKER_CLIENT_H
