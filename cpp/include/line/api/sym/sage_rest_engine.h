/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SYM_SAGE_REST_ENGINE_H
#define LINE_API_SYM_SAGE_REST_ENGINE_H

/**
 * @file
 * @ingroup api_sym
 * SymEngine backed by the line-sage-rest service.
 *
 * Port of jline.api.sym.SageRestEngine. The service is SageMath behind the JSON
 * protocol in io/sage/server.py. Every request is a single POST carrying the whole
 * problem, so nothing is bind-mounted and the client works against a container,
 * a remote host or a hand-started server alike.
 *
 * NUMERIC COEFFICIENTS ARE SENT AS DECIMAL STRINGS and read server side as
 * exact rationals, which is what keeps the solve exact: a double coerced by the
 * CAS would carry the binary rational nearest the decimal instead, and the
 * difference survives all the way into the printed normal form. The same rule
 * applies to the assignment eval() sends.
 *
 * The server enforces the timeout too, so a runaway symbolic solve is killed
 * there rather than merely abandoned here.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <vector>

#include "json.hpp"
#include "line/api/sym/sym_engine.h"
#include "line/util/error.h"
#include "line/util/http.h"

namespace line {
namespace sym {

namespace detail {

using Json = nlohmann::json;

/**
 * Canary for SageRestEngine::isUsable: the weighted-average softmin form, the
 * smallest expression observed to kill a worker whose FLINT wants BMI2/ADX on
 * a CPU that has neither. Its argument is a 17-digit decimal so the exact
 * rational is multi-limb, which is what reaches the offending routine.
 */
inline const char* canary_expr() {
    return "(x*exp(-x) + exp(-1))/(exp(-x) + exp(-1))";
}

/** @return the canary's argument */
inline const char* canary_arg() { return "0.68999999999999995"; }

/** The canary's value. */
const double CANARY_VALUE = 0.82116556904906557;

/** @return the usability verdicts, keyed by base URL */
inline std::map<std::string, bool>& usable_cache() {
    static std::map<std::string, bool> cache;
    return cache;
}

/** @return the mutex guarding usable_cache */
inline std::mutex& usable_mutex() {
    static std::mutex m;
    return m;
}

/**
 * The shortest decimal literal that round-trips to v, i.e. what Java's
 * Double.toString sends. Twin of line::reg::shortest_decimal, repeated here so
 * that api/sym stays independent of the numeric core it never otherwise uses.
 */
inline std::string decimal_string(double v) {
    char buf[64];
    for (int prec = 1; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof(buf), "%.*g", prec, v);
        if (std::strtod(buf, nullptr) == v) return std::string(buf);
    }
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

/** Reads a string array field, empty when the field is absent or null. */
inline std::vector<std::string> to_string_list(const Json& obj, const std::string& field) {
    std::vector<std::string> out;
    if (!obj.contains(field) || obj[field].is_null()) return out;
    for (const Json& el : obj[field]) out.push_back(el.get<std::string>());
    return out;
}

/** Reads a string field, or the fallback when it is absent or null. */
inline std::string opt_string(const Json& obj, const std::string& field,
                              const std::string& fallback) {
    if (obj.contains(field) && !obj[field].is_null()) return obj[field].get<std::string>();
    return fallback;
}

/** Encodes a square expression matrix, refusing a ragged one. */
inline Json to_json_matrix(const std::vector<std::vector<std::string>>& Q) {
    if (Q.empty()) throw InputError("SageRestEngine: Q must not be empty");
    Json rows = Json::array();
    for (std::size_t i = 0; i < Q.size(); ++i) {
        if (Q[i].size() != Q.size())
            throw InputError("SageRestEngine: Q must be square, row " + std::to_string(i) +
                             " has " + std::to_string(Q[i].size()) + " entries but Q has " +
                             std::to_string(Q.size()) + " rows");
        Json row = Json::array();
        for (std::size_t j = 0; j < Q[i].size(); ++j)
            row.push_back(Q[i][j].empty() ? std::string("0") : Q[i][j]);
        rows.push_back(row);
    }
    return rows;
}

/**
 * Encodes a string list, dropping empty entries. An empty symbol marks an event
 * with no positive rate, as symbolicGeneratorResult.symbols does in the JAR; it
 * contributes nothing and must not reach the server as "".
 */
inline Json to_json_array(const std::vector<std::string>& items) {
    Json arr = Json::array();
    for (std::size_t i = 0; i < items.size(); ++i)
        if (!items[i].empty()) arr.push_back(items[i]);
    return arr;
}

/**
 * A possibly NON-SQUARE block, e.g. a weight block.
 *
 * to_json_matrix is for the generator and insists on a square matrix; a weight
 * block is (number of measures) by (number of states) and is rarely square. The
 * column count is still checked, because a short row would otherwise be a
 * silently wrong measure rather than an error.
 */
inline Json to_json_rows(const std::vector<std::vector<std::string>>& rows, std::size_t width) {
    Json out = Json::array();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].size() != width)
            throw InputError("SageRestEngine: weight row " + std::to_string(i) + " has " +
                             std::to_string(rows[i].size()) + " entries but the generator has " +
                             std::to_string(width) + " states");
        Json row = Json::array();
        for (std::size_t j = 0; j < rows[i].size(); ++j)
            row.push_back(rows[i][j].empty() ? std::string("0") : rows[i][j]);
        out.push_back(row);
    }
    return out;
}

/**
 * A vector whose LENGTH IS PART OF THE CONTRACT, with an empty entry meaning
 * zero rather than absent.
 *
 * Deliberately not to_json_array: that one DROPS empty entries, which is right
 * for a symbol list (an event with no positive rate contributes no symbol) and
 * catastrophic for s0 or alpha, where dropping an entry shifts every state
 * after it and the server would see a shorter vector than the matrix.
 */
inline Json to_json_vector(const std::vector<std::string>& items) {
    Json arr = Json::array();
    for (std::size_t i = 0; i < items.size(); ++i)
        arr.push_back(items[i].empty() ? std::string("0") : items[i]);
    return arr;
}

inline std::vector<CtmcMeasure> to_measure_list(const Json& obj, const std::string& field) {
    std::vector<CtmcMeasure> out;
    if (!obj.contains(field) || obj[field].is_null()) return out;
    for (const Json& el : obj[field]) {
        CtmcMeasure m;
        m.name = opt_string(el, "name", "");
        m.expr = opt_string(el, "expr", "");
        m.num = opt_string(el, "num", "");
        m.den = opt_string(el, "den", "");
        m.reason = opt_string(el, "reason", "");
        out.push_back(m);
    }
    return out;
}

inline std::vector<int> to_int_list(const Json& obj, const std::string& field) {
    std::vector<int> out;
    if (!obj.contains(field) || obj[field].is_null()) return out;
    for (const Json& el : obj[field]) out.push_back(el.get<int>());
    return out;
}

}  // namespace detail

/** Client of the line-sage-rest service. */
class SageRestEngine : public SymEngine {
public:
    /** Default per-request timeout, in seconds. */
    static constexpr int DEFAULT_TIMEOUT_SECONDS = 300;

    /**
     * @param baseUrl base URL of the service, e.g. "http://localhost:8080"
     */
    explicit SageRestEngine(const std::string& baseUrl)
        : timeoutSeconds_(DEFAULT_TIMEOUT_SECONDS) {
        std::string u = baseUrl;
        const std::size_t b = u.find_first_not_of(" \t\r\n");
        const std::size_t e = u.find_last_not_of(" \t\r\n");
        u = b == std::string::npos ? std::string() : u.substr(b, e - b + 1);
        while (!u.empty() && u[u.size() - 1] == '/') u.erase(u.size() - 1);
        if (u.empty()) throw InputError("SageRestEngine: baseUrl must not be empty");
        baseUrl_ = u;
    }

    /**
     * Sets the per-request timeout.
     *
     * @param seconds timeout in seconds; not positive disables it
     * @return this engine
     */
    SageRestEngine& setTimeoutSeconds(int seconds) {
        timeoutSeconds_ = seconds;
        return *this;
    }

    /** @return the per-request timeout in seconds */
    int getTimeoutSeconds() const { return timeoutSeconds_; }

    /** @return the base URL this engine posts to */
    const std::string& getBaseUrl() const { return baseUrl_; }

    std::string name() const override { return "sage"; }

    bool isAvailable() const override {
        try {
            const detail::Json health = get("/api/v1/health", 5000);
            return detail::opt_string(health, "status", "") == "ok";
        } catch (const Error&) {
            return false;
        }
    }

    /**
     * Checks that the service can actually EVALUATE, not merely that it
     * answers.
     *
     * The line-sage-rest image ships a FLINT built for CPUs that have BMI2 and
     * ADX. On an older host the first multi-limb exact operation raises
     * SIGILL, the worker dies mid-request and the call returns no bytes at
     * all; /api/v1/health is pure Python and keeps answering, so it cannot see
     * this. The canary is the weighted-average softmin form, which is what the
     * fluid export actually sends, and is the smallest expression observed to
     * trigger it. Verdicts are cached per URL, so this costs one small request
     * the first time a service is considered and nothing after. See
     * _kb/11-conventions-and-gotchas.md.
     *
     * @return true if the service returned the canary's value
     */
    bool isUsable() const {
        {
            std::lock_guard<std::mutex> guard(detail::usable_mutex());
            std::map<std::string, bool>& cache = detail::usable_cache();
            const std::map<std::string, bool>::const_iterator it = cache.find(baseUrl_);
            if (it != cache.end()) return it->second;
        }
        bool ok = false;
        try {
            detail::Json request;
            request["exprs"] = detail::Json::array({std::string(detail::canary_expr())});
            detail::Json values = detail::Json::object();
            values["x"] = std::string(detail::canary_arg());
            request["values"] = values;
            request["timeout_s"] = 30;
            const detail::Json response =
                parse(http::post_json(baseUrl_ + "/api/v1/eval", request.dump(), 60000));
            checkStatus("/api/v1/eval", response);
            if (response.contains("values") && response["values"].is_array() &&
                response["values"].size() == 1 && !response["values"][0].is_null()) {
                const double v = response["values"][0].get<double>();
                ok = std::fabs(v - detail::CANARY_VALUE) < 1e-9;
            }
        } catch (const Error&) {
            // A dead worker closes the connection without a reply, which
            // surfaces as a transport error rather than a service one. Either
            // way the backend cannot serve us.
            ok = false;
        }
        if (!ok) {
            std::cerr << "[LINE] Ignoring symbolic backend at " << baseUrl_
                      << ": it did not return the usability canary. On a CPU without "
                      << "BMI2/ADX the image's FLINT raises SIGILL mid-request." << std::endl;
        }
        {
            std::lock_guard<std::mutex> guard(detail::usable_mutex());
            detail::usable_cache()[baseUrl_] = ok;
        }
        return ok;
    }

    /**
     * Reads the service identity, used to tell a line-sage-rest server apart
     * from another line-*-rest service on the same conventional port.
     *
     * @return the /api/v1/info document
     */
    detail::Json info() const { return get("/api/v1/info", 5000); }

    /**
     * Names the routes this client calls that the service does not serve.
     *
     * A SERVICE CAN BE HEALTHY AND STILL TOO OLD, and that case used to reach
     * the caller as a request failure rather than as a resolution failure. The
     * roster ran a line-sage-rest image predating the ctmc/measures and
     * ctmc/passage routes on 2026-09-11: it answered /api/v1/health, passed the
     * `isUsable` arithmetic canary, and then failed every call to those routes
     * with `unknown endpoint /api/v1/ctmc/measures` -- an error that reads as a
     * defect in the solver rather than as an out-of-date container.
     * /api/v1/info lists what the running server actually routes, so asking it
     * turns that into a named refusal at resolution time.
     *
     * @return the missing routes; empty when the service serves them all, and
     *         empty when it is too old to report its route list at all, which
     *         leaves the request itself to fail as before
     */
    std::vector<std::string> missingRoutes() const {
        static const char* kRequired[] = {
            "/api/v1/ctmc/solve", "/api/v1/ctmc/sensitivity", "/api/v1/ctmc/measures",
            "/api/v1/ctmc/passage", "/api/v1/simplify", "/api/v1/diff", "/api/v1/eval",
            "/api/v1/fluid/odes"};
        std::vector<std::string> missing;
        const detail::Json doc = info();
        if (!doc.contains("endpoints") || !doc["endpoints"].is_array()) return missing;
        const detail::Json& served = doc["endpoints"];
        for (std::size_t i = 0; i < sizeof(kRequired) / sizeof(kRequired[0]); ++i) {
            bool found = false;
            for (detail::Json::const_iterator it = served.begin(); it != served.end(); ++it) {
                if (it->is_string() && it->get<std::string>() == kRequired[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) missing.push_back(kRequired[i]);
        }
        return missing;
    }

    CtmcSolution solveCTMC(const std::vector<std::vector<std::string>>& Q,
                           const std::vector<std::string>& symbols) override {
        detail::Json request;
        request["Q"] = detail::to_json_matrix(Q);
        request["symbols"] = detail::to_json_array(symbols);
        request["normalize"] = true;
        const detail::Json response = post("/api/v1/ctmc/solve", request);

        CtmcSolution sol;
        sol.pi = detail::to_string_list(response, "pi");
        sol.num = detail::to_string_list(response, "num");
        sol.den = detail::opt_string(response, "den", "1");
        sol.nConnComp = response.contains("nConnComp") ? response["nConnComp"].get<int>() : 1;
        if (response.contains("connComp") && !response["connComp"].is_null())
            for (const detail::Json& el : response["connComp"]) sol.connComp.push_back(el.get<int>());
        return sol;
    }

    CtmcMeasures ctmcMeasures(const std::vector<std::vector<std::string>>& Q,
                              const std::vector<std::string>& symbols,
                              const std::vector<std::vector<std::string>>& weights,
                              const std::vector<std::string>& names,
                              const std::vector<RatioSpec>& ratios) override {
        detail::Json request;
        request["Q"] = detail::to_json_matrix(Q);
        request["symbols"] = detail::to_json_array(symbols);
        request["weights"] = detail::to_json_rows(weights, Q.size());
        request["normalize"] = true;
        if (!names.empty()) request["names"] = detail::to_json_vector(names);
        if (!ratios.empty()) {
            detail::Json arr = detail::Json::array();
            for (std::size_t i = 0; i < ratios.size(); ++i) {
                detail::Json o = detail::Json::object();
                o["name"] = ratios[i].name;
                o["num"] = ratios[i].num;
                o["den"] = ratios[i].den;
                arr.push_back(o);
            }
            request["ratios"] = arr;
        }
        const detail::Json response = post("/api/v1/ctmc/measures", request);

        CtmcMeasures out;
        out.pi = detail::to_string_list(response, "pi");
        out.num = detail::to_string_list(response, "num");
        out.den = detail::opt_string(response, "den", "1");
        out.nConnComp = response.contains("nConnComp") ? response["nConnComp"].get<int>() : 1;
        out.connComp = detail::to_int_list(response, "connComp");
        out.measures = detail::to_measure_list(response, "measures");
        out.ratios = detail::to_measure_list(response, "ratios");
        return out;
    }

    Passage ctmcPassage(const std::vector<std::vector<std::string>>& S,
                        const std::vector<std::string>& s0,
                        const std::vector<std::string>& alpha, const std::string& atom,
                        const std::vector<std::string>& symbols, const std::string& svar,
                        const std::vector<std::string>& want, int nmax) override {
        detail::Json request;
        request["S"] = detail::to_json_matrix(S);
        request["s0"] = detail::to_json_vector(s0);
        request["alpha"] = detail::to_json_vector(alpha);
        request["atom"] = atom.empty() ? std::string("0") : atom;
        request["symbols"] = detail::to_json_array(symbols);
        request["svar"] = svar.empty() ? std::string("s") : svar;
        request["want"] = detail::to_json_vector(want);
        request["nmax"] = nmax < 1 ? 1 : nmax;
        const detail::Json response = post("/api/v1/ctmc/passage", request);

        Passage out;
        out.lst = detail::opt_string(response, "lst", "");
        out.lstNum = detail::opt_string(response, "lst_num", "");
        out.lstDen = detail::opt_string(response, "lst_den", "");
        out.lstAll = detail::to_string_list(response, "lstall");
        out.moments = detail::to_string_list(response, "moments");
        if (response.contains("momall") && !response["momall"].is_null())
            for (const detail::Json& row : response["momall"]) {
                std::vector<std::string> r;
                for (const detail::Json& el : row) r.push_back(el.get<std::string>());
                out.momAll.push_back(r);
            }
        out.unreachable = detail::to_int_list(response, "unreachable");
        return out;
    }

    SymSensitivity ctmcSensitivity(const std::vector<std::vector<std::string>>& Q,
                                   const std::vector<std::string>& symbols,
                                   const std::string& theta,
                                   const std::vector<std::string>& reward) override {
        detail::Json request;
        request["Q"] = detail::to_json_matrix(Q);
        request["symbols"] = detail::to_json_array(symbols);
        request["theta"] = theta;
        if (!reward.empty()) request["reward"] = detail::to_json_array(reward);
        const detail::Json response = post("/api/v1/ctmc/sensitivity", request);

        SymSensitivity s;
        s.pi = detail::to_string_list(response, "pi");
        s.dpi = detail::to_string_list(response, "dpi");
        s.Er = detail::opt_string(response, "Er", "");
        s.S = detail::opt_string(response, "S", "");
        s.SS = detail::opt_string(response, "SS", "");
        s.hasReward = !reward.empty();
        return s;
    }

    std::vector<std::string> simplify(const std::vector<std::string>& exprs,
                                      const std::string& form) override {
        detail::Json request;
        request["exprs"] = detail::to_json_array(exprs);
        request["form"] = form.empty() ? std::string("cancel") : form;
        return detail::to_string_list(post("/api/v1/simplify", request), "results");
    }

    std::vector<std::string> diff(const std::vector<std::string>& exprs,
                                  const std::string& variable, int order) override {
        detail::Json request;
        request["exprs"] = detail::to_json_array(exprs);
        request["var"] = variable;
        request["order"] = order;
        return detail::to_string_list(post("/api/v1/diff", request), "results");
    }

    std::vector<double> eval(const std::vector<std::string>& exprs,
                             const std::map<std::string, double>& assignment) override {
        detail::Json request;
        request["exprs"] = detail::to_json_array(exprs);
        detail::Json values = detail::Json::object();
        for (std::map<std::string, double>::const_iterator it = assignment.begin();
             it != assignment.end(); ++it) {
            // sent as text so the server reads the decimal exactly, see the
            // header comment
            values[it->first] = detail::decimal_string(it->second);
        }
        request["values"] = values;
        const detail::Json response = post("/api/v1/eval", request);

        std::vector<double> out;
        if (!response.contains("values") || response["values"].is_null()) return out;
        for (const detail::Json& el : response["values"])
            out.push_back(el.is_null() ? std::nan("") : el.get<double>());
        return out;
    }

    FluidODEs fluidODEs(const std::vector<std::string>& rhs, const std::vector<std::string>& vars,
                        const std::vector<std::string>& want) override {
        detail::Json request;
        request["rhs"] = detail::to_json_array(rhs);
        request["vars"] = detail::to_json_array(vars);
        request["want"] = detail::to_json_array(want);
        const detail::Json response = post("/api/v1/fluid/odes", request);

        FluidODEs odes;
        if (response.contains("jacobian") && !response["jacobian"].is_null()) {
            odes.hasJacobian = true;
            for (const detail::Json& row : response["jacobian"]) {
                std::vector<std::string> r;
                for (const detail::Json& el : row) r.push_back(el.get<std::string>());
                odes.jacobian.push_back(r);
            }
        }
        if (response.contains("latex") && !response["latex"].is_null()) {
            odes.hasLatex = true;
            odes.latex = detail::to_string_list(response, "latex");
        }
        if (response.contains("equilibria") && !response["equilibria"].is_null()) {
            odes.hasEquilibria = true;
            for (const detail::Json& sol : response["equilibria"]) {
                std::map<std::string, std::string> m;
                for (detail::Json::const_iterator it = sol.begin(); it != sol.end(); ++it)
                    m[it.key()] = it.value().get<std::string>();
                odes.equilibria.push_back(m);
            }
        }
        return odes;
    }

private:
    detail::Json post(const std::string& path, detail::Json request) const {
        const int millis = timeoutSeconds_ > 0 ? timeoutSeconds_ * 1000 : 0;
        if (timeoutSeconds_ > 0) request["timeout_s"] = timeoutSeconds_;
        const detail::Json response = parse(
            http::post_json(baseUrl_ + path, request.dump(), millis));
        checkStatus(path, response);
        return response;
    }

    detail::Json get(const std::string& path, int millis) const {
        return parse(http::get(baseUrl_ + path, millis));
    }

    static detail::Json parse(const http::Response& response) {
        if (response.body.empty())
            throw SymEngineError("line-sage-rest returned HTTP " +
                                 std::to_string(response.status) + " with no body");
        detail::Json parsed = detail::Json::parse(response.body, nullptr, false);
        if (parsed.is_discarded() || !parsed.is_object())
            throw SymEngineError("line-sage-rest returned a non-JSON body: " + response.body);
        return parsed;
    }

    static void checkStatus(const std::string& path, const detail::Json& response) {
        const std::string status = detail::opt_string(response, "status", "");
        if (status == "ok") return;
        throw SymEngineError("line-sage-rest " + path + " failed [" +
                             detail::opt_string(response, "code", "error") +
                             "]: " + detail::opt_string(response, "message", "unspecified error"));
    }

    std::string baseUrl_;
    int timeoutSeconds_;
};

}  // namespace sym
}  // namespace line

#endif  // LINE_API_SYM_SAGE_REST_ENGINE_H
