/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_NETWORK_READER_H
#define LINE_IO_NETWORK_READER_H

/**
 * Reader for the LINE `model.json` interchange (a Network model) into a
 * `qn::Network<T>` built through the programmatic builder.
 *
 * The wire format is the one `linemodel_save.m`, the Python writer and the JAR
 * `LineModelIO` emit: a top-level `{format, version, model}` envelope whose
 * `model.type == "Network"` carries `nodes`, `classes` and `routing`. This
 * reader feeds the SAME builder + finalize the C++ programmatic API uses, so a
 * model that reaches C++ this way is indistinguishable from one authored in
 * code -- exactly the contract `lqn_builder.h` documents for its own reader.
 *
 * Scope: the queueing-network subset SolverMVA analyses. Distributions are
 * honoured on the moments the analyzers read (mean, SCV and, for QNA/polling,
 * the family the wire names). A construct outside this subset -- a
 * LayeredNetwork / Workflow / Environment model, an unsupported node kind, or a
 * distribution family whose moments this reader cannot reconstruct exactly --
 * is REFUSED BY NAME rather than silently degraded, matching the "a featset
 * name is a claim" rule the rest of the port follows.
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "json.hpp"
#include "line/lang/lang_types.h"
#include "line/lang/prior.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace detail {

using json = nlohmann::json;

/**
 * A model.json number that may arrive as the wire's INFINITY SPELLING.
 *
 * JSON has no infinity literal, so `linemodel_save` writes any infinite scalar
 * as the string `"Infinity"` / `"-Infinity"` -- a rule it applies to EVERY
 * numeric field, not to a named few. A plain `value("k", def)` therefore throws
 * `json.exception.type_error.302 type must be number, but is string` the moment
 * a model carries one, and the throw names no key: `spn_basic_open` (a
 * transition mode with `"numServers": "Infinity"`) failed identically under
 * every solver, which reads as a solver defect rather than as a reader gap.
 * The JAR (`LineModelIO`) and Python (`_servers_from_json`) readers both accept
 * the string form; this is the C++ twin.
 *
 * `nan` is also accepted, since the writer emits `null` for NaN and a reader
 * asked for a number from `null` would throw the same way.
 *
 * WHICH FIELDS NEED THIS, and why it is not every numeric read. The encoder's
 * rule is generic, but `linemodel_save` GUARDS the fields whose natural value
 * is infinite: `servers`, `buffer` and `classCap` are each emitted only when
 * `isfinite`, and `population` is written only for the closed and self-looping
 * classes, which are finite by construction. So the generic spelling is only
 * reachable on an UNGUARDED field. Transition-mode `numServers` is the proven
 * one, and all three writers agree it is the sentinel there. The remaining uses
 * below are the fields a non-MATLAB writer could still send infinite. Do not
 * take this as licence to convert the other `get<double>()` reads: a
 * distribution parameter that arrives as a string is a malformed model, and
 * should keep throwing.
 */
inline double num_from_json(const json& v) {
    if (v.is_null()) return std::numeric_limits<double>::quiet_NaN();
    if (v.is_string()) {
        const std::string s = v.get<std::string>();
        if (s == "Infinity" || s == "inf" || s == "Inf")
            return std::numeric_limits<double>::infinity();
        if (s == "-Infinity" || s == "-inf" || s == "-Inf")
            return -std::numeric_limits<double>::infinity();
        if (s == "NaN" || s == "nan") return std::numeric_limits<double>::quiet_NaN();
        // Anything else is a malformed number, not a zero: report it rather
        // than let atof turn a typo into a silently wrong model.
        const char* begin = s.c_str();
        char* end = nullptr;
        const double parsed = std::strtod(begin, &end);
        if (end == begin || *end != '\0')
            throw InputError("network_reader: expected a number, got the string '" + s + "'");
        return parsed;
    }
    if (v.is_boolean()) return v.get<bool>() ? 1.0 : 0.0;
    // Not a string: defer to nlohmann, so an object or an array still throws
    // the type error it always did.
    return v.get<double>();
}

/** `num_from_json` for an optional key, with the reader's default. */
inline double num_value(const json& obj, const char* key, double def) {
    return obj.contains(key) ? num_from_json(obj.at(key)) : def;
}

/**
 * The `dest` of a fork override record: a node NAME, or the empty string for
 * "every link this class takes", which the builder spells as node 0.
 *
 * An unknown name is an error and not a silent 0: dropping the destination
 * would spread an override meant for one link over all of them, which is a
 * different model that still solves.
 */
inline std::size_t fork_dest_index(const json& ov,
                                   const std::map<std::string, std::size_t>& node_idx,
                                   const std::string& fork_name) {
    if (!ov.contains("dest")) return 0;
    const std::string d = ov.at("dest").get<std::string>();
    if (d.empty()) return 0;
    const std::map<std::string, std::size_t>::const_iterator it = node_idx.find(d);
    if (it == node_idx.end())
        throw InputError("network_reader: the Fork node '" + fork_name +
                         "' declares an override towards '" + d + "', which is not a node");
    return it->second;
}

/** The `scheduling` attribute of a model.json Queue node. */
inline lang::SchedStrategy sched_from_json(const std::string& s) {
    using S = lang::SchedStrategy;
    if (s == "INF" || s == "inf") return S::INF;
    if (s == "FCFS" || s == "fcfs") return S::FCFS;
    if (s == "PS" || s == "ps") return S::PS;
    if (s == "LCFS" || s == "lcfs") return S::LCFS;
    if (s == "LCFSPR" || s == "lcfspr") return S::LCFSPR;
    if (s == "SIRO" || s == "siro") return S::SIRO;
    if (s == "HOL" || s == "hol") return S::HOL;
    if (s == "DPS" || s == "dps") return S::DPS;
    if (s == "GPS" || s == "gps") return S::GPS;
    if (s == "SEPT" || s == "sept") return S::SEPT;
    if (s == "LEPT" || s == "lept") return S::LEPT;
    if (s == "SJF" || s == "sjf") return S::SJF;
    if (s == "LJF" || s == "ljf") return S::LJF;
    if (s == "SRPT" || s == "srpt") return S::SRPT;
    if (s == "LPS" || s == "lps") return S::LPS;
    if (s == "POLLING" || s == "polling") return S::POLLING;
    // The preemptive, priority and pass-and-swap families. They were absent
    // here while `SchedStrategy` and the state layer carried them, so a model
    // the port can represent was refused at its own front door -- the wire
    // spells them exactly as the enum does, and `sched_to_text` is the list to
    // keep this one in step with.
    if (s == "FCFSPR" || s == "fcfspr") return S::FCFSPR;
    if (s == "FCFSPI" || s == "fcfspi") return S::FCFSPI;
    if (s == "FCFSPRIO" || s == "fcfsprio") return S::HOL;  // MATLAB's alias for HOL
    if (s == "FCFSPRPRIO" || s == "fcfsprprio") return S::FCFSPRPRIO;
    if (s == "FCFSPIPRIO" || s == "fcfspiprio") return S::FCFSPIPRIO;
    if (s == "LCFSPI" || s == "lcfspi") return S::LCFSPI;
    if (s == "LCFSPRIO" || s == "lcfsprio") return S::LCFSPRIO;
    if (s == "LCFSPRPRIO" || s == "lcfsprprio") return S::LCFSPRPRIO;
    if (s == "LCFSPIPRIO" || s == "lcfspiprio") return S::LCFSPIPRIO;
    if (s == "PSPRIO" || s == "psprio") return S::PSPRIO;
    if (s == "DPSPRIO" || s == "dpsprio") return S::DPSPRIO;
    if (s == "GPSPRIO" || s == "gpsprio") return S::GPSPRIO;
    if (s == "SRPT" || s == "srpt") return S::SRPT;
    if (s == "SRPTPRIO" || s == "srptprio") return S::SRPTPRIO;
    if (s == "PSJF" || s == "psjf") return S::PSJF;
    if (s == "FB" || s == "fb") return S::FB;
    if (s == "LRPT" || s == "lrpt") return S::LRPT;
    if (s == "SETF" || s == "setf") return S::SETF;
    if (s == "FSP" || s == "fsp") return S::FSP;
    if (s == "EDD" || s == "edd") return S::EDD;
    if (s == "EDF" || s == "edf") return S::EDF;
    if (s == "PAS" || s == "pas") return S::PAS;
    if (s == "OI" || s == "oi") return S::OI;
    if (s == "REF" || s == "ref") return S::REF;
    if (s == "EXT" || s == "ext") return S::EXT;
    throw UnsupportedError("network_reader: unsupported scheduling discipline '" + s + "'");
}

/** The `dropRule` string of a Queue node, matching the linemodel_io map. An
 * unrecognised value resolves to WAITQ, as MATLAB's str_to_droprule does. */
inline lang::DropStrategy drop_from_json(const std::string& s) {
    using D = lang::DropStrategy;
    if (s == "drop") return D::DROP;
    if (s == "waitingQueue") return D::WAITQ;
    if (s == "blockingAfterService") return D::BAS;
    if (s == "retrial") return D::RETRIAL;
    if (s == "retrialWithLimit") return D::RETRIAL_WITH_LIMIT;
    return D::WAITQ;
}

inline lang::ReplacementStrategy replacement_from_json(const std::string& s) {
    using R = lang::ReplacementStrategy;
    if (s == "RR" || s == "rr") return R::RR;
    if (s == "FIFO" || s == "fifo") return R::FIFO;
    if (s == "SFIFO" || s == "sfifo") return R::SFIFO;
    if (s == "LRU" || s == "lru") return R::LRU;
    if (s == "HLRU" || s == "hlru") return R::HLRU;
    if (s == "CLIMB" || s == "climb") return R::CLIMB;
    if (s == "QLRU" || s == "qlru") return R::QLRU;
    throw UnsupportedError("network_reader: unsupported cache replacement strategy '" + s + "'");
}

inline lang::PollingType polling_from_json(const std::string& s) {
    using P = lang::PollingType;
    if (s == "GATED" || s == "gated") return P::GATED;
    if (s == "EXHAUSTIVE" || s == "exhaustive") return P::EXHAUSTIVE;
    if (s == "KLIMITED" || s == "klimited" || s == "K-LIMITED") return P::KLIMITED;
    if (s == "DECREMENTING" || s == "decrementing") return P::DECREMENTING;
    throw UnsupportedError("network_reader: unsupported polling type '" + s + "'");
}

/**
 * A two-phase hyper-exponential matched to a mean and SCV (SCV >= 1), by the
 * balanced-means convention MATLAB `HyperExp.fitMeanAndSCV` and JMT use. Kept
 * local so the reader carries the exact same moments as the reference fit
 * rather than an approximation of them.
 */
template <class T>
lang::Distrib<T> hyperexp_fit_mean_scv(double mean, double scv) {
    if (!(scv >= 1.0))
        throw InputError("network_reader: HyperExp fitMeanAndSCV needs SCV >= 1, got " +
                         std::to_string(scv));
    const double p = 0.5 * (1.0 + std::sqrt((scv - 1.0) / (scv + 1.0)));
    const double mu1 = 2.0 * p / mean;
    const double mu2 = 2.0 * (1.0 - p) / mean;
    return lang::Distrib<T>::hyperexp(num_traits<T>::from_double(p),
                                      num_traits<T>::from_double(mu1),
                                      num_traits<T>::from_double(mu2));
}

/**
 * A Cache field, looked up in the FLAT spelling on the node first and then in
 * the nested `cache` object.
 *
 * The writers emit both (`linemodel_save.m:425-435`) so that the MATLAB and the
 * JAR readers each find what they look for; the two are mirrored from one
 * source and therefore agree. A free function rather than a lambda because the
 * result feeds `.get<...>()`, and a lambda declared inside a template makes
 * that call dependent.
 */
inline bool has_cache_key(const json& nd, const json& cj, const char* key) {
    return nd.contains(key) || cj.contains(key);
}
inline const json& cache_key(const json& nd, const json& cj, const char* key) {
    return nd.contains(key) ? nd.at(key) : cj.at(key);
}

/** The `routingStrategies` name of a dispatcher. */
inline lang::RoutingStrategy routing_from_json(const std::string& s) {
    typedef lang::RoutingStrategy R;
    if (s == "PROB") return R::PROB;
    if (s == "RAND") return R::RAND;
    if (s == "RROBIN") return R::RROBIN;
    if (s == "WRROBIN") return R::WRROBIN;
    if (s == "JSQ") return R::JSQ;
    if (s == "SQ" || s == "KCHOICES") return R::SQ;
    // The declaration itself rides in the `stateDepRouting` block; this only
    // names the entry row's strategy, which the block then supersedes.
    if (s == "SDR") return R::SDR;
    if (s == "FIRING") return R::FIRING;
    if (s == "DISABLED") return R::DISABLED;
    throw UnsupportedError("network_reader: unsupported routing strategy '" + s + "'");
}

/** The `patience` block's `impatienceType`. */
inline lang::ImpatienceType impatience_from_json(const std::string& s) {
    typedef lang::ImpatienceType I;
    if (s == "RENEGING") return I::RENEGING;
    if (s == "BALKING") return I::BALKING;
    if (s == "RETRIAL") return I::RETRIAL;
    throw UnsupportedError("network_reader: unsupported impatience type '" + s + "'");
}

/** The `balking` block's `strategy`. */
inline lang::BalkingStrategy balking_from_json(const std::string& s) {
    typedef lang::BalkingStrategy B;
    if (s == "QUEUE_LENGTH") return B::QUEUE_LENGTH;
    if (s == "EXPECTED_WAIT") return B::EXPECTED_WAIT;
    if (s == "COMBINED") return B::COMBINED;
    throw UnsupportedError("network_reader: unsupported balking strategy '" + s + "'");
}

/** The `heteroSchedPolicy` of a station with several server pools. */
inline lang::HeteroSchedPolicy hetero_from_json(const std::string& s) {
    typedef lang::HeteroSchedPolicy H;
    if (s == "ORDER") return H::ORDER;
    if (s == "ALIS") return H::ALIS;
    if (s == "ALFS") return H::ALFS;
    if (s == "FAIRNESS") return H::FAIRNESS;
    if (s == "FSF") return H::FSF;
    if (s == "RAIS") return H::RAIS;
    throw UnsupportedError("network_reader: unsupported heterogeneous scheduling policy '" + s + "'");
}

/** The `departureDiscipline` of a queueing Place. */
inline lang::DepartureDiscipline departure_from_json(const std::string& s) {
    if (s == "NORMAL") return lang::DepartureDiscipline::NORMAL;
    if (s == "FIFO") return lang::DepartureDiscipline::FIFO;
    throw UnsupportedError("network_reader: unsupported departure discipline '" + s + "'");
}

/** Read a wire field that may be a JSON array or, for length one, a bare scalar. */
template <class T>
std::vector<T> num_vec_from_json(const json& v) {
    std::vector<T> out;
    // num_from_json, not get<double>(): a state row is one of the UNGUARDED
    // fields of the comment above it. An open model's Source holds an infinite
    // population, so `initialState`/`stateSpace` carry "Infinity" the moment the
    // writer emits a state for every stateful node.
    if (v.is_array())
        for (const json& x : v) out.push_back(num_traits<T>::from_double(num_from_json(x)));
    else
        out.push_back(num_traits<T>::from_double(num_from_json(v)));
    return out;
}

/** A dense matrix written as an array of row arrays. */
template <class T>
Matrix<T> mat_from_json(const json& v) {
    if (!v.is_array())
        throw InputError("network_reader: expected an array of rows");
    Matrix<T> M(v.size(), v.empty() ? 0 : v.at(0).size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        const json& row = v.at(i);
        for (std::size_t j = 0; j < row.size(); ++j)
            M(i, j) = num_traits<T>::from_double(num_from_json(row.at(j)));
    }
    return M;
}

/** A list of dense matrices, the form the marked and batch families are written in. */
template <class T>
std::vector<Matrix<T> > mat_list_from_json(const json& v) {
    std::vector<Matrix<T> > out;
    for (const json& m : v) out.push_back(mat_from_json<T>(m));
    return out;
}

/** Reconstruct a distribution from a `fit` block (moments), by declared family. */
template <class T>
lang::Distrib<T> dist_from_fit(const std::string& type, const json& fit) {
    const std::string method = fit.at("method").get<std::string>();
    if (method == "fitMean" || method == "fitMeanAndOrder") {
        const double mean = fit.at("mean").get<double>();
        if (type == "Det") return lang::Distrib<T>::det(num_traits<T>::from_double(mean));
        if (type == "Erlang") {
            const long order = fit.contains("order") ? fit.at("order").get<long>() : 1L;
            const T rate = num_traits<T>::from_double(double(order) / mean);
            return lang::Distrib<T>::erlang(rate, std::size_t(order));
        }
        // fitMean on any other family is honoured on the mean alone -> Exp.
        return lang::Distrib<T>::exp_mean(num_traits<T>::from_double(mean));
    }
    if (method == "fitMeanAndSCV") {
        const double mean = fit.at("mean").get<double>();
        const double scv = fit.at("scv").get<double>();
        if (type == "Erlang")
            return lang::Distrib<T>::erlang_fit(num_traits<T>::from_double(mean),
                                                num_traits<T>::from_double(scv));
        if (type == "HyperExp") return hyperexp_fit_mean_scv<T>(mean, scv);
        if (std::fabs(scv - 1.0) < 1e-12)
            return lang::Distrib<T>::exp_mean(num_traits<T>::from_double(mean));
        throw UnsupportedError("network_reader: fitMeanAndSCV for family '" + type +
                               "' is not reconstructed; use explicit params");
    }
    throw UnsupportedError("network_reader: unsupported fit method '" + method + "' for '" + type +
                           "'");
}

/** Reconstruct a distribution from a model.json distribution record. */
template <class T>
lang::Distrib<T> dist_from_json(const json& obj) {
    const std::string type = obj.at("type").get<std::string>();
    if (type == "Immediate") return lang::Distrib<T>::immediate();
    if (type == "Disabled") return lang::Distrib<T>::disabled_dist();

    // A Prior arrives in one of its TWO forms, tagged by `kind` and defaulting
    // to the discrete one when the tag is absent, as every model.json written
    // before the continuous form existed is. The discrete form is the
    // alternative set plus the weights. The CONTINUOUS form is the parameter
    // density plus the distribution its factory BUILDS, with the parameter left
    // in the slots the factory fills (`linemodel_save.m:prior_factory2json`):
    // the handle itself cannot cross JSON, but substituting theta into those
    // slots reconstructs the same alternative at any node count, so
    // `options.samples` still means what it means in MATLAB. See
    // _kb/09-ldes-and-cache.md for the encoding and its one restriction.
    if (type == "Prior") {
        const std::string kind =
            obj.contains("kind") ? obj.at("kind").get<std::string>() : std::string("discrete");
        if (kind != "discrete" && kind != "continuous")
            throw InputError("network_reader: a Prior's 'kind' is 'discrete' or 'continuous', got '" +
                             kind + "'");
        if (kind == "continuous") {
            if (!obj.contains("paramDist") || !obj.contains("factory"))
                throw InputError(
                    "network_reader: a continuous Prior carries 'paramDist' and 'factory' (a "
                    "template distribution plus the parameter slots it fills)");
            const json& fac = obj.at("factory");
            if (!fac.contains("template") || !fac.contains("slots"))
                throw InputError(
                    "network_reader: a continuous Prior's 'factory' carries 'template' and "
                    "'slots'");
            const json tmpl = fac.at("template");
            std::vector<std::string> slots;
            for (const json& s : fac.at("slots")) slots.push_back(s.get<std::string>());
            if (slots.empty())
                throw InputError(
                    "network_reader: a continuous Prior's factory names no parameter slot, so the "
                    "parameter would not reach the distribution it builds");
            if (!tmpl.contains("params"))
                throw InputError(
                    "network_reader: a continuous Prior's factory template carries its parameters "
                    "as 'params'; a fitted or phase-type representation has no named slot");
            for (std::size_t l = 0; l < slots.size(); ++l)
                if (!tmpl.at("params").contains(slots[l]))
                    throw InputError("network_reader: a continuous Prior's factory names '" +
                                     slots[l] + "' as a parameter slot, but its template has no "
                                                "such parameter");
            // THETA CROSSES AS A DOUBLE, which costs nothing it was not already
            // costing: `dist_quantile` reaches it by bisection to FineTol, so the
            // parameter is an approximation of the stratum median well above
            // double resolution before this substitution sees it.
            const std::function<lang::Distrib<T>(const T&)> factory =
                [tmpl, slots](const T& theta) {
                    json j = tmpl;
                    for (std::size_t l = 0; l < slots.size(); ++l)
                        j["params"][slots[l]] = num_traits<T>::to_double(theta);
                    return dist_from_json<T>(j);
                };
            return lang::prior_continuous<T>(dist_from_json<T>(obj.at("paramDist")), factory);
        }
        if (!obj.contains("distributions") || !obj.contains("probabilities"))
            throw InputError(
                "network_reader: a discrete Prior carries 'distributions' and 'probabilities'");
        std::vector<lang::Distrib<T> > alts;
        for (const json& a : obj.at("distributions")) alts.push_back(dist_from_json<T>(a));
        std::vector<T> probs;
        for (const json& p : obj.at("probabilities"))
            probs.push_back(num_traits<T>::from_double(p.get<double>()));
        return lang::prior_discrete<T>(alts, probs);
    }

    // The time-inhomogeneous families. The wire form is the reference's
    // (`linemodel_save.m:2275-2302`, `linemodel_io.py:64-83`): a breakpoint
    // vector plus one matrix (or one row, or one scalar) per segment, and a
    // `cyclic` flag. NHPP carries `rates`, MAPt `D0`/`D1`, PHt `alpha`/`S`.
    if (type == "NHPP" || type == "MAPt" || type == "PHt") {
        const json& p = obj.at("params");
        std::vector<T> bp;
        for (const json& b : p.at("breakpoints"))
            bp.push_back(num_traits<T>::from_double(b.get<double>()));
        const bool cyc = p.contains("cyclic") && p.at("cyclic").get<bool>();
        auto read_mats = [&](const char* key) {
            std::vector<Matrix<T> > out;
            for (const json& seg : p.at(key)) {
                const std::vector<std::vector<double> > rows =
                    seg.get<std::vector<std::vector<double> > >();
                Matrix<T> M(rows.size(), rows.empty() ? 0 : rows[0].size());
                for (std::size_t a = 0; a < rows.size(); ++a)
                    for (std::size_t b = 0; b < rows[a].size(); ++b)
                        M(a, b) = num_traits<T>::from_double(rows[a][b]);
                out.push_back(M);
            }
            return out;
        };
        if (type == "NHPP") {
            std::vector<T> rates;
            for (const json& r : p.at("rates"))
                rates.push_back(num_traits<T>::from_double(r.get<double>()));
            return lang::Distrib<T>::nhpp(bp, rates, cyc);
        }
        if (type == "MAPt")
            return lang::Distrib<T>::mapt(bp, read_mats("D0"), read_mats("D1"), cyc);
        std::vector<std::vector<T> > alphas;
        for (const json& a : p.at("alpha")) {
            const std::vector<double> row = a.get<std::vector<double> >();
            std::vector<T> av;
            for (double v : row) av.push_back(num_traits<T>::from_double(v));
            alphas.push_back(av);
        }
        return lang::Distrib<T>::pht(bp, alphas, read_mats("S"), cyc);
    }

    // A Replayer names a trace FILE, and the writers add the APH fit of its
    // first three moments beside it precisely because that path may not resolve
    // on another machine. Reading the trace back gives the same distribution
    // MATLAB had; the fit is the documented fallback, and the stored mean the
    // last resort. This branch precedes the `ph` one below, which would
    // otherwise turn every Replayer into a bare PH and drop the trace.
    if (type == "Replayer" || type == "Trace") {
        if (obj.contains("params") && obj.at("params").contains("fileName")) {
            const std::string path = obj.at("params").at("fileName").get<std::string>();
            std::ifstream tr(path.c_str());
            if (tr) {
                std::vector<T> samples;
                double v = 0.0;
                while (tr >> v) samples.push_back(num_traits<T>::from_double(v));
                if (!samples.empty()) {
                    lang::Distrib<T> d = lang::Distrib<T>::replayer(samples);
                    d.trace_file = path;  // only the exporters read it back
                    return d;
                }
            }
        }
        if (obj.contains("ph")) {
            const json& ph = obj.at("ph");
            const std::vector<double> a = ph.at("alpha").get<std::vector<double> >();
            const std::vector<std::vector<double> > rows =
                ph.at("T").get<std::vector<std::vector<double> > >();
            std::vector<T> alpha;
            for (double x : a) alpha.push_back(num_traits<T>::from_double(x));
            Matrix<T> A(rows.size(), rows.empty() ? 0 : rows[0].size());
            for (std::size_t i = 0; i < rows.size(); ++i)
                for (std::size_t j = 0; j < rows[i].size(); ++j)
                    A(i, j) = num_traits<T>::from_double(rows[i][j]);
            return lang::Distrib<T>::phase_type(alpha, A, true);
        }
        if (obj.contains("params") && obj.at("params").contains("mean"))
            return lang::Distrib<T>::exp_mean(
                num_traits<T>::from_double(obj.at("params").at("mean").get<double>()));
        throw InputError(
            "network_reader: a Replayer carries neither a readable trace file, nor the APH fit "
            "the writers add beside it, nor a mean; there is nothing to reconstruct");
    }

    // Markovian arrival families carry an explicit (D0, D1). A generic MAP/MMAP
    // arrives as a `map` object; an MMPP2 as its four rate params.
    if (type == "MMPP2") {
        const json& p = obj.at("params");
        const double l0 = p.at("lambda0").get<double>(), l1 = p.at("lambda1").get<double>();
        const double s0 = p.at("sigma0").get<double>(), s1 = p.at("sigma1").get<double>();
        Matrix<T> D0(2, 2), D1(2, 2);
        D1(0, 0) = num_traits<T>::from_double(l0);
        D1(1, 1) = num_traits<T>::from_double(l1);
        D0(0, 0) = num_traits<T>::from_double(-(l0 + s0));
        D0(0, 1) = num_traits<T>::from_double(s0);
        D0(1, 0) = num_traits<T>::from_double(s1);
        D0(1, 1) = num_traits<T>::from_double(-(l1 + s1));
        return lang::Distrib<T>::map_dist(D0, D1, lang::ProcessType::MMPP2);
    }
    if (type == "MAP" || type == "MMPP") {
        if (!obj.contains("map"))
            throw InputError("network_reader: " + type + " carries no 'map' (D0, D1) object");
        const json& mp = obj.at("map");
        return lang::Distrib<T>::map_dist(mat_from_json<T>(mp.at("D0")),
                                          mat_from_json<T>(mp.at("D1")),
                                          lang::ProcessType::MAP);
    }
    // A MARKED MAP carries one arrival block per mark in its own `mmap` object;
    // the aggregate D1 the unmarked consumers read is their sum, rebuilt on load
    // exactly as the reference readers do.
    if (type == "MMAP" || type == "MarkedMAP") {
        if (!obj.contains("mmap"))
            throw InputError("network_reader: " + type + " carries no 'mmap' (D0, D1k) object");
        const json& mp = obj.at("mmap");
        return lang::Distrib<T>::mmap(mat_from_json<T>(mp.at("D0")),
                                      mat_list_from_json<T>(mp.at("D1k")));
    }
    // A BMAP writes the WHOLE block list D0, D1, ..., Dk, where Dj carries an
    // arrival of batch size j; a MarkedMMPP writes the same list plus the mark
    // count K, and lowers to the MMAP process type as MATLAB's fromText does.
    if (type == "BMAP" || type == "MarkedMMPP") {
        const json& p = obj.at("params");
        const std::vector<Matrix<T> > D = mat_list_from_json<T>(p.at("D"));
        if (type == "BMAP") return lang::Distrib<T>::bmap(D);
        if (D.size() < 2)
            throw InputError("network_reader: a MarkedMMPP carries D0 and at least one marked block");
        return lang::Distrib<T>::mmap(D[0], std::vector<Matrix<T> >(D.begin() + 1, D.end()));
    }
    if (type == "DMAP") {
        const json& p = obj.at("params");
        return lang::Distrib<T>::dmap(mat_from_json<T>(p.at("D0")), mat_from_json<T>(p.at("D1")));
    }
    if (type == "RAP") {
        const json& p = obj.at("params");
        return lang::Distrib<T>::rap(mat_from_json<T>(p.at("H0")), mat_from_json<T>(p.at("H1")));
    }
    // ME and CME share one process type: the concentrated representation is a
    // subclass carrying no distinct tag, exactly as MATLAB `fromText` maps it.
    if (type == "ME" || type == "CME") {
        const json& p = obj.at("params");
        const std::vector<double> a = p.at("alpha").get<std::vector<double> >();
        std::vector<T> alpha;
        for (double v : a) alpha.push_back(num_traits<T>::from_double(v));
        return lang::Distrib<T>::me(alpha, mat_from_json<T>(p.at("A")));
    }

    // A phase-type family carries its representation as `ph: {alpha, T}`, which
    // is the third form the writer emits alongside `params` and `fit`. Without
    // this branch a PH/APH/Coxian written with an explicit representation --
    // exactly what `Coxian.fit`, `APH.fitMeanAndSCV` and the MAP-to-PH paths
    // produce -- was refused as "has neither params nor a fit block", which is a
    // statement about what the reader looked for and not about what the writer
    // sent. The representation is complete on the wire; only the branch was
    // missing.
    if (obj.contains("ph")) {
        const json& ph = obj.at("ph");
        const std::vector<double> a = ph.at("alpha").get<std::vector<double> >();
        const std::vector<std::vector<double> > rows =
            ph.at("T").get<std::vector<std::vector<double> > >();
        if (a.empty() || rows.size() != a.size())
            throw InputError("network_reader: distribution '" + type +
                             "' has a 'ph' block whose alpha and T disagree in order");
        std::vector<T> alpha;
        alpha.reserve(a.size());
        for (double v : a) alpha.push_back(num_traits<T>::from_double(v));
        Matrix<T> A(rows.size(), rows.empty() ? 0 : rows[0].size());
        for (std::size_t i = 0; i < rows.size(); ++i)
            for (std::size_t j = 0; j < rows[i].size(); ++j)
                A(i, j) = num_traits<T>::from_double(rows[i][j]);
        // APH and PH share a representation and differ only in the type tag,
        // so the wire's own name decides it rather than a structural test.
        return lang::Distrib<T>::phase_type(alpha, A, type == "APH");
    }

    // Explicit parameters take precedence over a fit block, mirroring the
    // reference readers.
    if (!obj.contains("params") && obj.contains("fit"))
        return dist_from_fit<T>(type, obj.at("fit"));
    if (!obj.contains("params"))
        throw InputError("network_reader: distribution '" + type +
                         "' has neither params nor a fit block");
    const json& p = obj.at("params");
    if (type == "Exp") {
        const double lam = p.contains("lambda") ? p.at("lambda").get<double>()
                                                : p.at("rate").get<double>();
        return lang::Distrib<T>::exp_rate(num_traits<T>::from_double(lam));
    }
    if (type == "Det") return lang::Distrib<T>::det(num_traits<T>::from_double(p.at("value").get<double>()));
    if (type == "Erlang") {
        const double lam = p.at("lambda").get<double>();
        const long k = p.at("k").get<long>();
        return lang::Distrib<T>::erlang(num_traits<T>::from_double(lam), std::size_t(k));
    }
    if (type == "HyperExp") {
        const std::vector<T> pv = num_vec_from_json<T>(p.at("p"));
        const std::vector<T> lv = num_vec_from_json<T>(p.at("lambda"));
        // The two-branch entry point carries MATLAB's getParam order, which a
        // parameter dump compares against; wider ones take the general form.
        if (pv.size() == 2 && lv.size() == 2)
            return lang::Distrib<T>::hyperexp(pv[0], lv[0], lv[1]);
        return lang::Distrib<T>::hyperexp_n(pv, lv);
    }
    if (type == "Coxian")
        return lang::Distrib<T>::coxian(num_vec_from_json<T>(p.at("mu")),
                                        num_vec_from_json<T>(p.at("phi")));
    if (type == "Cox2")
        return lang::Distrib<T>::cox2(num_traits<T>::from_double(p.at("mu1").get<double>()),
                                      num_traits<T>::from_double(p.at("mu2").get<double>()),
                                      num_traits<T>::from_double(p.at("phi1").get<double>()));
    if (type == "Uniform")
        return lang::Distrib<T>::uniform(num_traits<T>::from_double(p.at("a").get<double>()),
                                         num_traits<T>::from_double(p.at("b").get<double>()));
    // Gamma(alpha = shape, beta = SCALE), the pairing the JAR loader fixes.
    if (type == "Gamma")
        return lang::Distrib<T>::gamma_dist(num_traits<T>::from_double(p.at("alpha").get<double>()),
                                            num_traits<T>::from_double(p.at("beta").get<double>()));
    if (type == "Lognormal")
        return lang::Distrib<T>::lognormal(num_traits<T>::from_double(p.at("mu").get<double>()),
                                           num_traits<T>::from_double(p.at("sigma").get<double>()));
    // `Normal` shares the (mu, sigma) key pair with `Lognormal` above, and is
    // the one family here that is never a service process: it crosses the wire
    // only as the parameter density of a continuous Prior. Refusing it made
    // every such Prior a model MATLAB could write and this reader could not
    // read, which is a UQ study that stops at the interchange.
    if (type == "Normal")
        return lang::Distrib<T>::normal(num_traits<T>::from_double(p.at("mu").get<double>()),
                                        num_traits<T>::from_double(p.at("sigma").get<double>()));
    // Pareto(alpha = shape, scale) and Weibull(alpha = SCALE, beta = SHAPE):
    // the two families spell the wire key `alpha` for opposite roles, which is
    // the reference's convention and the one trap in this table.
    if (type == "Pareto")
        return lang::Distrib<T>::pareto(num_traits<T>::from_double(p.at("alpha").get<double>()),
                                        num_traits<T>::from_double(p.at("scale").get<double>()));
    if (type == "Weibull")
        return lang::Distrib<T>::weibull(num_traits<T>::from_double(p.at("alpha").get<double>()),
                                         num_traits<T>::from_double(p.at("beta").get<double>()));
    if (type == "DiscreteUniform")
        return lang::Distrib<T>::discrete_uniform(
            num_traits<T>::from_double(p.at("min").get<double>()),
            num_traits<T>::from_double(p.at("max").get<double>()));
    if (type == "Bernoulli")
        return lang::Distrib<T>::bernoulli(num_traits<T>::from_double(p.at("p").get<double>()));
    if (type == "Binomial")
        return lang::Distrib<T>::binomial(num_traits<T>::from_double(p.at("n").get<double>()),
                                          num_traits<T>::from_double(p.at("p").get<double>()));
    if (type == "Poisson")
        return lang::Distrib<T>::poisson(num_traits<T>::from_double(p.at("lambda").get<double>()));
    if (type == "Geometric")
        return lang::Distrib<T>::geometric(num_traits<T>::from_double(p.at("p").get<double>()));
    if (type == "Zipf")
        return lang::Distrib<T>::zipf(num_traits<T>::from_double(p.at("s").get<double>()),
                                      std::size_t(p.at("n").get<long>()));
    if (type == "DiscreteSampler") {
        const std::vector<T> pv = num_vec_from_json<T>(p.at("p"));
        const std::vector<T> xv =
            p.contains("x") ? num_vec_from_json<T>(p.at("x")) : std::vector<T>();
        return lang::Distrib<T>::discrete_sampler(pv, xv);
    }
    if (type == "EmpiricalCDF" || type == "EmpiricalCdf")
        return lang::Distrib<T>::empirical_cdf(num_vec_from_json<T>(p.at("x")),
                                               num_vec_from_json<T>(p.at("F")));

    // THE MOMENT-ONLY FALLBACK the writers emit for a family with no JSON
    // representation of its own: `{type: <name>, params: {mean, scv?}}`, written
    // with a warning at save time. It is honoured on the moments it carries --
    // an exponential from a mean alone, and otherwise the acyclic phase-type
    // matching the two -- rather than refused, because refusing would make a
    // model unreadable that the reference itself declares readable. The TYPE
    // TAG IS NOT KEPT: the reconstruction is a different law that happens to
    // share two moments, and claiming the original name for it would make
    // `sn.procid` lie about what every solver is actually integrating.
    if (p.contains("mean") && !p.contains("scv"))
        return lang::Distrib<T>::exp_mean(num_traits<T>::from_double(p.at("mean").get<double>()));
    if (p.contains("mean") && p.contains("scv")) {
        const double mean = p.at("mean").get<double>();
        const double scv = p.at("scv").get<double>();
        if (std::fabs(scv - 1.0) < 1e-12)
            return lang::Distrib<T>::exp_mean(num_traits<T>::from_double(mean));
        if (scv < 1.0)
            return lang::Distrib<T>::erlang_fit(num_traits<T>::from_double(mean),
                                                num_traits<T>::from_double(scv));
        return hyperexp_fit_mean_scv<T>(mean, scv);
    }
    throw UnsupportedError("network_reader: unsupported distribution family '" + type +
                           "' (params form)");
}

/**
 * The item-popularity pmf of a cache read class, `nodeparam.pread`.
 *
 * A popularity is a DISCRETE distribution over the item ranks, and the writers
 * emit whichever family the model used: `DiscreteSampler` carries the pmf
 * itself, `Zipf` carries only (s, n) and the pmf p_i = i^-s / H(s,n) is
 * rebuilt here. Reading only the first form left every Zipf-popularity cache
 * with an EMPTY read vector, which is a uniform cache by default rather than a
 * refusal.
 */
/**
 * The popularity LAW, recorded beside the pmf for the exporters.
 *
 * JMT's Cache section takes a parametric popularity and cannot be given a pmf,
 * so the exponent and the support size have to survive the read; every solver
 * still reads the pmf `pmf_from_json` returns.
 */
template <class T>
typename qn::CacheParam<T>::Popularity popularity_kind_from_json(const json& obj,
                                                                 std::size_t nitems) {
    typename qn::CacheParam<T>::Popularity k;
    const std::string type = obj.value("type", std::string());
    const json& p = obj.contains("params") ? obj.at("params") : obj;
    if (type == "Zipf") {
        k.type = lang::ProcessType::ZIPF;
        k.s = p.at("s").get<double>();
        k.n = p.contains("n") ? std::size_t(p.at("n").get<long>()) : nitems;
    } else if (type == "DiscreteSampler" || p.contains("p")) {
        k.type = lang::ProcessType::DISCRETESAMPLER;
        k.n = p.contains("p") ? p.at("p").size() : nitems;
    }
    return k;
}

template <class T>
std::vector<T> pmf_from_json(const json& obj, std::size_t nitems) {
    const std::string type = obj.value("type", std::string());
    const json& p = obj.contains("params") ? obj.at("params") : obj;
    std::vector<T> out;
    // THE ARRAY TEST IS THE WHOLE GUARD. `params.p` is also the SCALAR success
    // probability of Bernoulli, Binomial and Geometric, and `p.contains("p")`
    // alone admitted those into this branch, where iterating a JSON scalar
    // yields one element and the family collapsed to a one-point pmf.
    if (type == "DiscreteSampler" || (p.contains("p") && p.at("p").is_array())) {
        for (const json& v : p.at("p")) out.push_back(num_traits<T>::from_double(v.get<double>()));
        return out;
    }
    if (type == "Zipf") {
        const double s = p.at("s").get<double>();
        const std::size_t n = p.contains("n") ? std::size_t(p.at("n").get<long>()) : nitems;
        double h = 0.0;
        for (std::size_t k = 1; k <= n; ++k) h += std::pow(double(k), -s);
        for (std::size_t k = 1; k <= n; ++k)
            out.push_back(num_traits<T>::from_double(std::pow(double(k), -s) / h));
        return out;
    }
    throw UnsupportedError(
        "network_reader: a cache popularity is written as '" + type +
        "', and the discrete families carrying an item pmf are DiscreteSampler and Zipf");
}

/**
 * How far an infinitely-supported batch law is materialized.
 *
 * `sn.signalremdist` is a VECTOR, so a Geometric or a Poisson has to stop
 * somewhere, and the consumer (`signal_batch_pmf`) lumps everything past the
 * last stored entry onto "remove the whole eligible population". That lumping
 * is the REFERENCE'S OWN rule only for mass beyond the population present, so
 * the vector must outrun any population a station can hold; 4096 does, for a
 * chain whose per-class cutoff is a few hundred at the very most.
 */
constexpr std::size_t REMOVAL_PMF_MAX_TERMS = 4096;

/**
 * The batch-size pmf of a G-network signal, INDEXED BY THE BATCH SIZE ITSELF.
 *
 * WHY THIS IS NOT `pmf_from_json`. That function returns a cache popularity,
 * whose entry i is the mass of RANK i+1; here entry b is the mass of the batch
 * size b, counting from a batch of ZERO. `sn.signalremdist` is documented on the
 * second convention (network_struct.h) and both consumers -- `signal_batch_pmf`
 * for the chain and the LDES draw -- read it that way, so a reader that returns
 * the first shifts every batch by one job.
 *
 * MATLAB EVALUATES THE LAW, IT DOES NOT TABULATE IT: `State.signalBatchPMF`
 * calls `dist.evalPMF(0:ntot)`, so the wire carries the FAMILY (`{type:
 * "Geometric", params: {p: 0.5}}`) and the pmf has to be rebuilt here. Reading
 * only `DiscreteSampler` sent every parametric family through a branch that
 * iterated the scalar `params.p` and produced the one-point pmf
 * P(B = 0) = p -- a signal that removes NOTHING with probability p and empties
 * the station otherwise. On the `test_batch_removal` G-network that read
 * utilization 0.58690 against the reference's 0.56422.
 */
template <class T>
std::vector<T> removal_pmf_from_json(const json& obj) {
    const std::string type = obj.value("type", std::string());
    const json& p = obj.contains("params") ? obj.at("params") : obj;
    std::vector<T> out;
    auto put = [&out](std::size_t k, double v) {
        if (out.size() <= k) out.resize(k + 1, num_traits<T>::from_int(0));
        out[k] = num_traits<T>::from_double(num_traits<T>::to_double(out[k]) + v);
    };
    if (type == "DiscreteSampler" || (p.contains("p") && p.at("p").is_array())) {
        // `x` IS THE SUPPORT AND NOT A LABEL: MATLAB's `DiscreteSampler(p, x)`
        // defaults it to 1..n and `evalPMF` looks a value up in it, so dropping
        // it would read the first mass as a batch of zero.
        const json& pv = p.at("p");
        const bool has_x = p.contains("x") && p.at("x").is_array();
        for (std::size_t i = 0; i < pv.size(); ++i) {
            const double xi =
                has_x ? p.at("x").at(i).get<double>() : static_cast<double>(i + 1);
            if (xi < 0) throw InputError("network_reader: a batch size cannot be negative");
            put(static_cast<std::size_t>(xi + 0.5), pv.at(i).get<double>());
        }
        return out;
    }
    if (type == "Bernoulli") {
        const double q = p.at("p").get<double>();
        put(0, 1.0 - q);
        put(1, q);
        return out;
    }
    if (type == "Binomial") {
        const double q = p.at("p").get<double>();
        const std::size_t n = static_cast<std::size_t>(p.at("n").get<double>() + 0.5);
        double term = std::pow(1.0 - q, static_cast<double>(n));
        for (std::size_t k = 0; k <= n; ++k) {
            put(k, term);
            if (k < n && q < 1.0)
                term *= (static_cast<double>(n - k) / static_cast<double>(k + 1)) * q / (1.0 - q);
        }
        return out;
    }
    if (type == "Poisson") {
        const double lam = p.at("lambda").get<double>();
        double term = std::exp(-lam), acc = 0.0;
        for (std::size_t k = 0; k < REMOVAL_PMF_MAX_TERMS; ++k) {
            put(k, term);
            acc += term;
            if (acc > 1.0 - 1e-15) break;
            term *= lam / static_cast<double>(k + 1);
        }
        return out;
    }
    if (type == "Geometric") {
        // MATLAB's Geometric(p) counts TRIALS TO THE FIRST SUCCESS, support
        // {1, 2, ...}: P(B = 0) is zero, so a signal always takes at least one.
        const double q = p.at("p").get<double>();
        if (q <= 0.0 || q > 1.0) throw InputError("network_reader: Geometric(p) needs 0 < p <= 1");
        put(0, 0.0);
        double term = q, acc = 0.0;
        for (std::size_t k = 1; k <= REMOVAL_PMF_MAX_TERMS; ++k) {
            put(k, term);
            acc += term;
            if (acc > 1.0 - 1e-15) break;
            term *= (1.0 - q);
        }
        return out;
    }
    if (type == "DiscreteUniform") {
        const long lo = static_cast<long>(p.at("min").get<double>());
        const long hi = static_cast<long>(p.at("max").get<double>());
        if (hi < lo || lo < 0)
            throw InputError("network_reader: DiscreteUniform(min,max) needs 0 <= min <= max");
        const double w = 1.0 / static_cast<double>(hi - lo + 1);
        for (long k = lo; k <= hi; ++k) put(static_cast<std::size_t>(k), w);
        return out;
    }
    if (type == "Det") {
        const double v = p.contains("t") ? p.at("t").get<double>() : p.at("mean").get<double>();
        if (v < 0) throw InputError("network_reader: a batch size cannot be negative");
        put(static_cast<std::size_t>(v + 0.5), 1.0);
        return out;
    }
    throw UnsupportedError(
        "network_reader: a signal removal law is written as '" + type +
        "', and the discrete families a batch size can be drawn from are DiscreteSampler, "
        "Bernoulli, Binomial, Poisson, Geometric, DiscreteUniform and Det");
}

/**
 * Rebuild a class- or joint-dependence handle from the box-lattice table the
 * writers emit.
 *
 * A function handle cannot cross JSON, so `linemodel_save.m` (`cd_scaling_table`),
 * the JAR `LineModelIO` and the Python writer all MATERIALIZE beta(n) / eta_i(n)
 * over 0 <= n(r) <= cutoffs(r), keyed by the comma-joined 0-based per-class
 * counts. This rebuilds the callable from that table, clamping the population to
 * the cutoffs so the scaling SATURATES beyond the tabulated range exactly as the
 * table intends, and returning all-ones for a composition absent from the table
 * so an unlisted state leaves the nominal rate unscaled. Twin of the Python
 * `_cd_table_to_callable`.
 */
template <class T>
lang::CdScaling<T> cd_scaling_from_json(const json& tbl, const std::vector<int>& cutoffs,
                                        std::size_t K) {
    std::map<std::string, std::vector<T> > table;
    for (auto it = tbl.begin(); it != tbl.end(); ++it)
        table[it.key()] = num_vec_from_json<T>(it.value());
    const std::vector<int> cut = cutoffs;
    return [table, cut, K](const std::vector<T>& n) {
        std::string key;
        for (std::size_t r = 0; r < K; ++r) {
            int v = r < n.size()
                        ? static_cast<int>(std::lround(num_traits<T>::to_double(n[r])))
                        : 0;
            if (v < 0) v = 0;
            if (r < cut.size() && v > cut[r]) v = cut[r];
            if (r) key += ',';
            key += std::to_string(v);
        }
        std::vector<T> out(K, num_traits<T>::from_int(1));
        typename std::map<std::string, std::vector<T> >::const_iterator it = table.find(key);
        if (it == table.end()) return out;
        for (std::size_t r = 0; r < K && r < it->second.size(); ++r) out[r] = it->second[r];
        return out;
    };
}

/**
 * Rebuild the mu(c) of an order-independent / pass-and-swap station from the
 * macrostate table the writers materialize (`oi_rate_table` in
 * `linemodel_save.m`, `LineModelIO.oiServiceRate` in the JAR).
 *
 * The handle takes an ORDERED microstate -- the list of 1-based class indices
 * in buffer order, which is what `set_pas` is defined over -- and the table is
 * keyed by the per-class COUNTS, because mu is order-independent by
 * construction. A composition beyond the tabulated cutoffs saturates at them,
 * matching how the writer chose those cutoffs (the closed populations, or 10
 * for an open class beyond which mu is constant); a composition the table does
 * not list returns zero, which is the writer's own encoding of a non-finite
 * rate and means the macrostate is unreachable.
 */
template <class T>
std::function<T(const std::vector<std::size_t>&)> oi_rate_from_json(const json& tbl,
                                                                    const std::vector<int>& cutoffs,
                                                                    std::size_t K) {
    std::map<std::string, T> table;
    for (auto it = tbl.begin(); it != tbl.end(); ++it)
        table[it.key()] = num_traits<T>::from_double(it.value().get<double>());
    const std::vector<int> cut = cutoffs;
    return [table, cut, K](const std::vector<std::size_t>& micro) {
        std::vector<int> cnt(K, 0);
        for (std::size_t j = 0; j < micro.size(); ++j)
            if (micro[j] >= 1 && micro[j] <= K) ++cnt[micro[j] - 1];
        std::string key;
        for (std::size_t r = 0; r < K; ++r) {
            int v = cnt[r];
            if (r < cut.size() && v > cut[r]) v = cut[r];
            if (r) key += ',';
            key += std::to_string(v);
        }
        typename std::map<std::string, T>::const_iterator it = table.find(key);
        return it == table.end() ? num_traits<T>::from_int(0) : it->second;
    };
}

/**
 * The declared peak the wire carries, or, for legacy JSON written before the
 * peak became mandatory, the peak DERIVED from the table.
 *
 * Deriving it is not a fabrication: the table is the whole lattice the reference
 * `cd_peak_scaling` would sweep, so the maximum over its rows and classes is
 * exactly that routine's answer. It skips the all-zero composition (the
 * reference's `tot > 0` guard) and any non-finite entry, which a handle may
 * legitimately return for an unreachable composition and which would otherwise
 * become the normalizer and zero every utilization at the station.
 */
template <class T>
std::vector<T> cd_peak_from_json(const json& blk, const json& tbl) {
    if (blk.contains("peak") && !blk.at("peak").empty())
        return num_vec_from_json<T>(blk.at("peak"));
    double bmax = 0;
    for (auto it = tbl.begin(); it != tbl.end(); ++it) {
        if (it.key().find_first_not_of("0,") == std::string::npos) continue;
        const std::vector<double> row = num_vec_from_json<double>(it.value());
        for (double x : row)
            if (std::isfinite(x) && x > bmax) bmax = x;
    }
    return std::vector<T>(1, num_traits<T>::from_double(bmax));
}

/**
 * Refuse any top-level or node-level key that carries model semantics this
 * reader does not consume.
 *
 * The whitelists are the keys the passes below actually read, plus the purely
 * descriptive ones. Anything else is a construct the C++ model layer either
 * cannot represent or does not yet parse, and it is NAMED rather than dropped.
 */
inline void reject_unconsumed_model_keys(const json& model) {
    static const char* kModelKeys[] = {"name",    "type",    "nodes",
                                       "classes", "routing", "format",
                                       "version", "rewards", "finiteCapacityRegions",
                                       "routingStrategies", "routingWeights", "routingParams",
                                       "logPath", "globalDependence", "stateDepRouting"};
    // EXACTLY the keys a branch below reads, and nothing else. The first draft
    // of this list also carried `arrival`, `classCap`, `accessGraph`,
    // `joinStrategy`, `joinQuorum`, `fanOut` and `swapGraph`,
    // none of which this reader consumes -- whitelisting them would have
    // re-admitted the very silent drop the gate exists to stop. When adding a
    // key here, grep that the parser actually reads it.
    static const char* kNodeKeys[] = {
        "name",       "type",        "scheduling",   "servers",     "service",
        "buffer",     "capacity",    "dropRule",     "schedParams", "pollingType",
        "pollingPar", "classSwitchMatrix",           "csMatrix",    "forkNode",
        "tasksPerLink",              "fanOutByDest", "fanOutDist",  "fanOutProb",
        "items",      "numItems",    "itemLevelCap", "popularity",  "replacementStrategy",
        "itemSizes",  "costCaps",    "accessProb",  "admissionProb", "accessGraph",
        "itemClass",
        "cache",      "initialState",
        "retrievalSystem",           "queues",       "hitClass",    "missClass",
        "immediateFeedback",         "loadDependence",
        "classDependence",           "jointDependence",
        "modes",      "classCap",    "departureDiscipline",
        "oiServiceRate",             "oiCutoffs",    "swapGraph",
        "arrivalBatch",              "markedClasses",
        "stateSpace", "statePrior",  "joinStrategy", "joinQuorum",
        "setupTime",  "delayOffTime","switchoverTimes",  "breakdown",
        "serverTypes","heteroSchedPolicy",     "serverParallelism",
        "balking",    "retrial",     "patience",     "orbitImpatience",
        "batchRejectProb",
        // Logger trace configuration; read into NodeDef::logger below.
        "fileName",   "filePath",    "startTime",    "loggerName",  "timestamp",
        "jobID",      "jobClass",    "timeSameClass","timeAnyClass"};
    // The class object. It had NO gate until every key below was found to be
    // read: `isReferenceClass`, `deadline`, `patience`, `spawnClass` and
    // `replySignalClass` were all being dropped in silence, each of them a
    // different wrong number rather than a diagnostic.
    static const char* kClassKeys[] = {
        "name",        "type",         "population",   "refNode",      "priority",
        "openOrClosed","signalType",   "targetClass",  "removalPolicy","removalDistribution",
        "isReferenceClass",            "deadline",     "patience",     "impatienceType",
        "spawnClass",  "replySignalClass",             "immediateFeedback"};
    auto known = [](const char* const* tab, std::size_t n, const std::string& k) {
        for (std::size_t i = 0; i < n; ++i)
            if (k == tab[i]) return true;
        return false;
    };
    const std::string why =
        "', which this reader does not implement. Refusing rather than dropping it: a "
        "constraint silently discarded here would make every solver return a confident "
        "answer for a different model";
    for (auto it = model.begin(); it != model.end(); ++it)
        if (!known(kModelKeys, sizeof(kModelKeys) / sizeof(*kModelKeys), it.key()))
            throw UnsupportedError("network_reader: the model carries '" + it.key() + why);
    if (model.contains("classes"))
        for (const json& cl : model.at("classes"))
            for (auto it = cl.begin(); it != cl.end(); ++it)
                if (!known(kClassKeys, sizeof(kClassKeys) / sizeof(*kClassKeys), it.key()))
                    throw UnsupportedError("network_reader: class '" +
                                           cl.value("name", std::string("?")) + "' carries '" +
                                           it.key() + why);
    if (!model.contains("nodes")) return;
    for (const json& nd : model.at("nodes"))
        for (auto it = nd.begin(); it != nd.end(); ++it)
            if (!known(kNodeKeys, sizeof(kNodeKeys) / sizeof(*kNodeKeys), it.key()))
                throw UnsupportedError("network_reader: node '" +
                                       nd.value("name", std::string("?")) + "' carries '" +
                                       it.key() + why);
}

}  // namespace detail

/**
 * Build a `qn::Network<T>` from a parsed model.json envelope.
 *
 * Two passes: classes and nodes are declared first (a class references its
 * reference node by name, a Join references its Fork by name, so every node
 * index must exist before the cross-references are wired), then service,
 * arrival and routing are applied.
 */
template <class T>
qn::Network<T> build_network_from_json(const detail::json& root) {
    using detail::json;
    const json& model = root.contains("model") ? root.at("model") : root;
    const std::string mtype = model.value("type", std::string("Network"));
    if (mtype != "Network") {
        // An Environment IS solved by this port, just not by this reader, so
        // the refusal names the arm that reads it rather than leaving the
        // caller to conclude the model is unsupported.
        if (mtype == "Environment")
            throw UnsupportedError(
                "network_reader: this is an Environment model (a network per stage plus the "
                "stage transitions); solve it with -s env, which reads it through "
                "environment_reader.h");
        throw UnsupportedError("network_reader: model type '" + mtype +
                               "' is not a Network; only Network models are solved by this path");
    }

    // WHAT THIS READER DOES NOT UNDERSTAND, IT REFUSES. A key carrying model
    // semantics that no branch below consumes would otherwise be SILENTLY
    // DROPPED, and every solver would then return an answer correct for the
    // model received and wrong for the model intended, with no diagnostic
    // anywhere. Not hypothetical: `fcr_mm1kdrop` exports a
    // `finiteCapacityRegions` block with globalMaxJobs 3 and a drop rule, and
    // without this gate the C++ MVA reported QLen 4 -- exactly rho/(1-rho) for
    // the UNBOUNDED M/M/1 -- against the exact M/M/1/K value 1.224932. Same
    // rule, and same reason, as the unknown-argument refusal on the --api
    // boundary: a key that silently takes its default is a wrong answer.
    detail::reject_unconsumed_model_keys(model);

    qn::Network<T> net(model.value("name", std::string("model")));
    net.set_log_path(model.value("logPath", std::string()));

    // A ClassSwitch node needs the class count at construction and a Join needs
    // its Fork to exist, while a closed class needs its reference node. That
    // cycle is broken by ordering the declarations, not by post-hoc setters:
    // every node a class or a Join can depend on is created first (pass 1a),
    // then the classes (1b), then the dependent nodes (1c). Routing is wired by
    // node NAME (pass 3), so this reordering never perturbs the parity table.
    const json& nodes = model.at("nodes");
    std::map<std::string, std::size_t> node_idx;
    std::vector<std::string> node_type(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i)
        node_type[i] = nodes[i].at("type").get<std::string>();

    // -- Pass 1a: nodes that carry no class/fork dependency. -----------------
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const json& nd = nodes[i];
        const std::string& type = node_type[i];
        // A Join is a STATION, so deferring its creation would shift every
        // station index after it and rotate the station rows of the result
        // document; it is declared here, in model order, and bound to its Fork
        // in pass 1c. ClassSwitch, Cache and Transition are plain nodes, so
        // their late creation perturbs no station index.
        if (type == "ClassSwitch" || type == "Cache" || type == "Transition") continue;
        const std::string name = nd.at("name").get<std::string>();
        std::size_t idx = 0;
        if (type == "Source") {
            idx = net.add_source(name);
        } else if (type == "Sink") {
            idx = net.add_sink(name);
        } else if (type == "Delay") {
            idx = net.add_delay(name);
        } else if (type == "Queue") {
            idx = net.add_queue(name,
                                detail::sched_from_json(nd.value("scheduling", std::string("FCFS"))));
        } else if (type == "Router") {
            idx = net.add_router(name);
        } else if (type == "Logger" || type == "LogTunnel") {
            idx = net.add_logger(name, nd.value("fileName", std::string()));
            qn::NodeDef::LoggerParam& lg = net.raw_struct().nodes[idx - 1].logger;
            if (nd.contains("filePath")) lg.file_path = nd.at("filePath").get<std::string>();
            lg.start_time = nd.value("startTime", lg.start_time);
            lg.logger_name = nd.value("loggerName", lg.logger_name);
            lg.timestamp = nd.value("timestamp", lg.timestamp);
            lg.job_id = nd.value("jobID", lg.job_id);
            lg.job_class = nd.value("jobClass", lg.job_class);
            lg.time_same_class = nd.value("timeSameClass", lg.time_same_class);
            lg.time_any_class = nd.value("timeAnyClass", lg.time_any_class);
        } else if (type == "Place") {
            idx = net.add_place(name);
        } else if (type == "Join") {
            idx = net.add_join_unbound(name);
        } else if (type == "Fork") {
            idx = net.add_fork(name, detail::num_value(nd, "tasksPerLink", 1.0));
        } else {
            throw UnsupportedError("network_reader: unsupported node type '" + type + "' at node '" +
                                   name + "'");
        }
        node_idx[name] = idx;
    }

    // -- Pass 1b: classes. ---------------------------------------------------
    const json& classes = model.at("classes");
    std::map<std::string, std::size_t> class_idx;
    // The signal classes, resolved after the whole class list exists: a
    // signal's TARGET is another class by name, which may follow it.
    std::vector<std::size_t> signal_classes;
    for (std::size_t r = 0; r < classes.size(); ++r) {
        const json& cl = classes[r];
        const std::string name = cl.at("name").get<std::string>();
        const std::string type = cl.at("type").get<std::string>();
        std::size_t idx = 0;
        if (type == "Open") {
            idx = net.add_open_class(name, cl.value("priority", 0));
        } else if (type == "Closed" || type == "SelfLooping") {
            // A `SelfLoopingClass` IS a closed class -- its jobs cycle at the
            // reference station, which is a property of the routing -- so it is
            // built as one and only tagged, exactly as `linemodel_load.m:125`
            // handles the two labels in one branch.
            const double pop = cl.at("population").get<double>();
            const std::string ref = cl.at("refNode").get<std::string>();
            auto it = node_idx.find(ref);
            if (it == node_idx.end())
                throw InputError("network_reader: class '" + name + "' references unknown node '" +
                                 ref + "'");
            idx = type == "SelfLooping"
                      ? net.add_self_looping_class(name, pop, it->second, cl.value("priority", 0))
                      : net.add_closed_class(name, pop, it->second, cl.value("priority", 0));
        } else if (type == "Signal") {
            // A G-network SIGNAL is an ordinary open or closed class that
            // REMOVES jobs instead of joining a queue, so it is created as its
            // underlying kind first and marked afterwards. `openOrClosed`
            // carries that kind; without it the signal would be built as a
            // closed class with no population.
            const std::string kind = cl.value("openOrClosed", std::string("Open"));
            if (kind == "Closed") {
                const std::string ref = cl.at("refNode").get<std::string>();
                auto it = node_idx.find(ref);
                if (it == node_idx.end())
                    throw InputError("network_reader: signal class '" + name +
                                     "' references unknown node '" + ref + "'");
                // A CLOSED SIGNAL HAS NO POPULATION OF ITS OWN, and the writers
                // therefore emit no `population` for it: `ClosedSignal.m:69`
                // passes 0 up to ClosedClass, its jobs arriving only by class
                // switch from the caller. Requiring the key made every REPLY
                // model unreadable here with a raw json out_of_range.
                idx = net.add_closed_class(name, cl.value("population", 0.0), it->second,
                                           cl.value("priority", 0));
            } else {
                idx = net.add_open_class(name, cl.value("priority", 0));
            }
            signal_classes.push_back(r);
        } else {
            throw UnsupportedError("network_reader: unsupported class type '" + type +
                                   "' for class '" + name + "'");
        }
        // The class-wide spelling of immediate feedback; the node-level map is
        // read in pass 2 and `sn.immfeed` is the OR of the two.
        if (cl.value("immediateFeedback", false)) net.set_class_immediate_feedback(idx);
        // `setReferenceClass`: which class of a chain `sn.refclass` names. It is
        // the denominator of every chain visit ratio in `sn_get_demands_chain`
        // and `sn_get_product_form_params`, so dropping it silently rescales
        // the demands of every multi-class chain that declares one.
        if (cl.value("isReferenceClass", false)) net.set_reference_class(idx);
        // `JobClass.deadline`, `sn.classdeadline`: EDD and EDF order by it.
        if (cl.contains("deadline")) {
            const double due = cl.at("deadline").get<double>();
            if (std::isfinite(due)) net.set_class_deadline(idx, due);
        }
        // The CLASS-WIDE patience, `Queue.getPatience`'s fallback. A node-scoped
        // `patience` read in pass 2 overrides it at the station that names it.
        if (cl.contains("patience")) {
            const json& pt = cl.at("patience");
            const lang::Distrib<T> pd = detail::dist_from_json<T>(pt);
            if (!pd.disabled)
                net.set_class_patience(
                    idx, pd,
                    cl.contains("impatienceType")
                        ? detail::impatience_from_json(cl.at("impatienceType").get<std::string>())
                        : lang::ImpatienceType::RENEGING);
        }
        class_idx[name] = idx;
    }
    // The class-to-class bindings, resolved once every class exists: either
    // side may be declared after the class that names it.
    for (std::size_t r = 0; r < classes.size(); ++r) {
        const json& cl = classes[r];
        const std::size_t idx = class_idx.at(cl.at("name").get<std::string>());
        if (cl.contains("spawnClass")) {
            const std::string sp = cl.at("spawnClass").get<std::string>();
            auto sit = class_idx.find(sp);
            if (sit == class_idx.end())
                throw InputError("network_reader: class '" + cl.at("name").get<std::string>() +
                                 "' spawns class '" + sp + "', which the model does not declare");
            net.set_class_spawn(idx, sit->second);
        }
        // Without this a REPLY signal class is INERT after a round trip:
        // nothing unblocks the servers waiting on it (`linemodel_save.m:797`).
        if (cl.contains("replySignalClass")) {
            const std::string rp = cl.at("replySignalClass").get<std::string>();
            auto rit = class_idx.find(rp);
            if (rit == class_idx.end())
                throw InputError("network_reader: class '" + cl.at("name").get<std::string>() +
                                 "' replies with class '" + rp +
                                 "', which the model does not declare");
            net.set_reply_signal_class(idx, rit->second);
        }
    }
    // Signals, now that every class the target may name exists.
    for (std::size_t si = 0; si < signal_classes.size(); ++si) {
        const json& cl = classes[signal_classes[si]];
        const std::string name = cl.at("name").get<std::string>();
        const std::string st = cl.value("signalType", std::string("negative"));
        lang::SignalType kind = lang::SignalType::NEGATIVE;
        if (st == "reply" || st == "REPLY") kind = lang::SignalType::REPLY;
        else if (st == "catastrophe" || st == "CATASTROPHE") kind = lang::SignalType::CATASTROPHE;
        else if (st != "negative" && st != "NEGATIVE")
            throw UnsupportedError("network_reader: signal class '" + name + "' is of type '" + st +
                                   "', and the kinds on the wire are negative, catastrophe and "
                                   "reply");
        lang::RemovalPolicy pol = lang::RemovalPolicy::RANDOM;
        const std::string rp = cl.value("removalPolicy", std::string("RANDOM"));
        if (rp == "FCFS" || rp == "fcfs") pol = lang::RemovalPolicy::FCFS;
        else if (rp == "LCFS" || rp == "lcfs") pol = lang::RemovalPolicy::LCFS;
        std::size_t target = 0;
        if (cl.contains("targetClass")) {
            auto tit = class_idx.find(cl.at("targetClass").get<std::string>());
            if (tit == class_idx.end())
                throw InputError("network_reader: signal class '" + name + "' targets class '" +
                                 cl.at("targetClass").get<std::string>() +
                                 "', which the model does not declare");
            target = tit->second;
        }
        std::vector<T> remdist;
        if (cl.contains("removalDistribution"))
            remdist = detail::removal_pmf_from_json<T>(cl.at("removalDistribution"));
        net.set_signal(class_idx.at(name), kind, pol, target, remdist);
    }
    const std::size_t K = class_idx.size();

    // -- Pass 1c: ClassSwitch (matrix now sizeable), Join (fork exists) and
    //    Cache (its hit/miss classes and popularity are class-indexed). --------
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const json& nd = nodes[i];
        const std::string& type = node_type[i];
        if (type != "ClassSwitch" && type != "Join" && type != "Cache" && type != "Fork") continue;
        const std::string name = nd.at("name").get<std::string>();
        std::size_t idx = 0;
        if (type == "Fork") {
            // VARIABLE FORKING LEVELS. Each list is an array of
            // {dest, class, ...} records, exactly as `linemodel_save.m:509-553`
            // writes them; `dest` is a node NAME and `class` a 1-based index.
            // The overrides are recorded here and replayed by `link()` in Pass 3,
            // because a per-destination override reads the routing.
            if (!nd.contains("fanOutByDest") && !nd.contains("fanOutDist") &&
                !nd.contains("fanOutProb"))
                continue;
            idx = node_idx.at(name);
            if (nd.contains("fanOutByDest")) {
                const json& ovs = nd.at("fanOutByDest");
                for (std::size_t e = 0; e < ovs.size(); ++e)
                    net.set_fork_tasks_per_link(idx, ovs[e].at("class").get<std::size_t>(),
                                                ovs[e].at("value").get<double>(),
                                                detail::fork_dest_index(ovs[e], node_idx, name));
            }
            if (nd.contains("fanOutDist")) {
                const json& ovs = nd.at("fanOutDist");
                for (std::size_t e = 0; e < ovs.size(); ++e) {
                    const std::vector<double> pv = ovs[e].at("p").get<std::vector<double> >();
                    const std::vector<double> xv = ovs[e].at("x").get<std::vector<double> >();
                    std::vector<T> p, x;
                    for (std::size_t q = 0; q < pv.size(); ++q)
                        p.push_back(num_traits<T>::from_double(pv[q]));
                    for (std::size_t q = 0; q < xv.size(); ++q)
                        x.push_back(num_traits<T>::from_double(xv[q]));
                    net.set_fork_tasks_per_link_dist(idx, ovs[e].at("class").get<std::size_t>(),
                                                     lang::Distrib<T>::discrete_sampler(p, x),
                                                     detail::fork_dest_index(ovs[e], node_idx, name));
                }
            }
            if (nd.contains("fanOutProb")) {
                const json& ovs = nd.at("fanOutProb");
                for (std::size_t e = 0; e < ovs.size(); ++e)
                    net.set_fork_branch_probability(
                        idx, ovs[e].at("class").get<std::size_t>(),
                        detail::fork_dest_index(ovs[e], node_idx, name),
                        ovs[e].at("value").get<double>());
            }
            continue;
        }
        if (type == "Cache") {
            // THE WRITERS EMIT THE CACHE TWICE, by design: a nested `cache`
            // object and the same fields flattened onto the node, so that both
            // the MATLAB and the JAR readers find what they look for
            // (`linemodel_save.m:425-435`). The flat spelling wins where both
            // are present -- they are mirrored from the nested one, so they
            // agree -- and the nested one is consulted for anything the flat
            // mirror does not carry.
            const json empty_obj = json::object();
            const json& cj = nd.contains("cache") ? nd.at("cache") : empty_obj;
            qn::CacheParam<T> cp;
            cp.nitems = detail::has_cache_key(nd, cj, "numItems")
                            ? detail::cache_key(nd, cj, "numItems").get<std::size_t>()
                            : cj.at("items").get<std::size_t>();
            cp.itemcap = detail::has_cache_key(nd, cj, "itemLevelCap")
                             ? detail::cache_key(nd, cj, "itemLevelCap").get<std::vector<int> >()
                             : cj.at("capacity").get<std::vector<int> >();
            // Per-item storage costs and per-list cost caps (ton21cache Sec. IX).
            // A scalar costCaps is the single cache-wide cap, replicated per list.
            if (detail::has_cache_key(nd, cj, "itemSizes"))
                cp.itemsize = detail::cache_key(nd, cj, "itemSizes").get<std::vector<int> >();
            if (detail::has_cache_key(nd, cj, "costCaps")) {
                const json& cc = detail::cache_key(nd, cj, "costCaps");
                if (cc.is_array()) {
                    cp.costcap = cc.get<std::vector<int> >();
                } else {
                    cp.costcapglobal = true;
                    cp.costcap.assign(cp.itemcap.size(), cc.get<int>());
                }
            }
            cp.replacestrat = detail::replacement_from_json(
                detail::has_cache_key(nd, cj, "replacementStrategy")
                    ? detail::cache_key(nd, cj, "replacementStrategy").get<std::string>()
                    : cj.value("replacement", std::string("RR")));
            cp.pread.assign(K, std::vector<T>());
            cp.hitclass.assign(K, 0);
            cp.missclass.assign(K, 0);
            cp.preadkind.assign(K, typename qn::CacheParam<T>::Popularity());
            cp.classitem.assign(K, 0);
            if (detail::has_cache_key(nd, cj, "popularity")) {
                const json& pop = detail::cache_key(nd, cj, "popularity");
                for (auto it = pop.begin(); it != pop.end(); ++it) {
                    // A class that does not read the cache is written with a
                    // `Disabled` popularity, and leaves `pread` empty
                    // (`linemodel_load.m:602`).
                    if (it.value().value("type", std::string()) == "Disabled") continue;
                    cp.pread[class_idx.at(it.key()) - 1] =
                        detail::pmf_from_json<T>(it.value(), cp.nitems);
                    cp.preadkind[class_idx.at(it.key()) - 1] =
                        detail::popularity_kind_from_json<T>(it.value(), cp.nitems);
                }
            }
            auto fill_switch = [&](const char* key, std::vector<std::size_t>& dst) {
                if (!detail::has_cache_key(nd, cj, key)) return;
                const json& blk = detail::cache_key(nd, cj, key);
                for (auto it = blk.begin(); it != blk.end(); ++it)
                    dst[class_idx.at(it.key()) - 1] = class_idx.at(it.value().get<std::string>());
            };
            fill_switch("hitClass", cp.hitclass);
            fill_switch("missClass", cp.missclass);
            // `itemClass`: for a cache network, the item each per-item class reads
            // (`Cache.setItemReadClasses`). Carried so a reader need not infer it
            // from a one-hot pread, which a genuine single-item popularity also has.
            if (detail::has_cache_key(nd, cj, "itemClass")) {
                const json& blk = detail::cache_key(nd, cj, "itemClass");
                for (auto it = blk.begin(); it != blk.end(); ++it)
                    cp.classitem[class_idx.at(it.key()) - 1] =
                        static_cast<std::size_t>(it.value().get<double>());
            }
            // `accessProb`: the per-(class, item) access graph, `sn.nodeparam{i}.accost`.
            // The wire form is a class-major array of item-major arrays of
            // (h+1)x(h+1) matrices, with an empty entry where the pair declares
            // none (`linemodel_save.m:402-416`). It is READ HERE and nowhere else,
            // and dropping it would solve a cache whose admission and promotion
            // are item-dependent as if they were the linear chain.
            // `admissionProb`: the q-LRU admission probability, `nodeparam.qlru`.
            // Every Cache the reference writes carries it, defaulting to 1, so
            // refusing it made every MATLAB-written cache model unreadable here.
            if (detail::has_cache_key(nd, cj, "admissionProb"))
                cp.qlru = num_traits<T>::from_double(detail::cache_key(nd, cj, "admissionProb").get<double>());
            if (detail::has_cache_key(nd, cj, "accessProb")) {
                const json& ap = detail::cache_key(nd, cj, "accessProb");
                cp.accost.clear();
                for (const json& per_class : ap) {
                    std::vector<Matrix<T> > row;
                    for (const json& g : per_class) {
                        if (g.is_null() || g.empty()) {
                            row.push_back(Matrix<T>());
                            continue;
                        }
                        row.push_back(detail::mat_from_json<T>(g));
                    }
                    cp.accost.push_back(row);
                }
            } else if (detail::has_cache_key(nd, cj, "accessGraph")) {
                // `accessGraph` is the SAME access cost shared by every class:
                // one (h+1)x(h+1) matrix per item, written when the model set
                // `node.graph` rather than the full per-class `accessProb`
                // (`linemodel_save.m:395-401`). Replicating it across the classes
                // here is what makes the two spellings the same model.
                const std::vector<Matrix<T> > shared =
                    detail::mat_list_from_json<T>(detail::cache_key(nd, cj, "accessGraph"));
                cp.accost.assign(K, shared);
            }
            // The initial cache contents: the state row `[class counts | list
            // contents | retrieval bitmap]` the reference dumps from the node.
            if (detail::has_cache_key(nd, cj, "initialState"))
                cp.initstate = detail::num_vec_from_json<T>(detail::cache_key(nd, cj, "initialState"));
            // Delayed-hit retrieval system: the retrieval classes are already in
            // the class list (created in pass 1b) with their service and routing,
            // so only the cache-side maps are recovered here -- no re-creation.
            // Read through the same nested-or-flat helper every other cache key
            // uses. All three writers put this one flat today, so a direct
            // `nd.contains` happens to work -- but it is the shape python got
            // wrong (it chose ONE source per node and lost whatever the other
            // held), and a reader that treats one key differently from its
            // siblings is the trap that made that possible.
            if (detail::has_cache_key(nd, cj, "retrievalSystem")) {
                const json& rs = detail::cache_key(nd, cj, "retrievalSystem");
                cp.retrieval_capacity = rs.value("capacity", 0);
                cp.retrieval_classes.assign(cp.nitems, std::vector<std::size_t>(K, 0));
                if (rs.contains("byClass")) {
                    for (auto rc = rs.at("byClass").begin(); rc != rs.at("byClass").end(); ++rc) {
                        const std::size_t rdcls = class_idx.at(rc.key());  // 1-based read class
                        std::vector<std::size_t> qnodes;
                        for (const auto& qn : rc.value().at("queues"))
                            qnodes.push_back(node_idx.at(qn.get<std::string>()));
                        cp.retrieval_queues[rdcls - 1] = qnodes;
                        for (auto it2 = rc.value().at("items").begin();
                             it2 != rc.value().at("items").end(); ++it2) {
                            const std::size_t item = std::stoul(it2.key());
                            if (item < cp.nitems)
                                cp.retrieval_classes[item][rdcls - 1] =
                                    class_idx.at(it2.value().get<std::string>());
                        }
                    }
                }
            }
            idx = net.add_cache(name, cp);
            node_idx[name] = idx;
            continue;
        }
        if (type == "Join") {
            std::size_t fork = 0;
            if (nd.contains("forkNode")) {
                auto it = node_idx.find(nd.at("forkNode").get<std::string>());
                if (it == node_idx.end())
                    throw InputError("network_reader: Join '" + name +
                                     "' references unknown fork '" +
                                     nd.at("forkNode").get<std::string>() + "'");
                fork = it->second;
            }
            // The station itself was declared in pass 1a, in model order, so
            // only the Fork it closes is recorded here.
            idx = node_idx.at(name);
            net.bind_join(idx, fork);
            // The declared join rule. PARTIAL fires on a quorum of siblings
            // rather than on all of them, so dropping it would make a partial
            // join wait for arrivals that never come.
            if (nd.contains("joinStrategy") || nd.contains("joinQuorum")) {
                const std::string js = nd.value("joinStrategy", std::string("STD"));
                if (js != "STD" && js != "PARTIAL")
                    throw UnsupportedError("network_reader: Join '" + name + "' declares strategy '" +
                                           js + "', and the rules on the wire are STD and PARTIAL");
                net.set_join_strategy(idx,
                                      js == "PARTIAL" ? lang::JoinStrategy::PARTIAL
                                                      : lang::JoinStrategy::STD,
                                      nd.value("joinQuorum", 0.0));
            }
        } else {
            // Rows a class-switch matrix omits leave that class unchanged, so the
            // matrix defaults to the identity before the listed rows overwrite it.
            // The builder's class indices are 1-based; the matrix is 0-indexed.
            // `csMatrix` is the JAR writer's spelling of the same block.
            Matrix<T> C(K, K);
            for (std::size_t d = 0; d < K; ++d) C(d, d) = num_traits<T>::from_int(1);
            if (nd.contains("classSwitchMatrix") || nd.contains("csMatrix")) {
                const json& csm =
                    nd.contains("classSwitchMatrix") ? nd.at("classSwitchMatrix") : nd.at("csMatrix");
                for (auto ri = csm.begin(); ri != csm.end(); ++ri) {
                    const std::size_t rr = class_idx.at(ri.key()) - 1;
                    for (std::size_t d = 0; d < K; ++d) C(rr, d) = num_traits<T>::from_int(0);
                    for (auto ci = ri.value().begin(); ci != ri.value().end(); ++ci)
                        C(rr, class_idx.at(ci.key()) - 1) =
                            num_traits<T>::from_double(ci.value().get<double>());
                }
            }
            idx = net.add_class_switch(name, C);
        }
        node_idx[name] = idx;
    }

    // -- Pass 1d: Transitions, whose arcs name other nodes. -------------------
    // A mode's enabling/inhibiting/firing arcs are written as (node, class,
    // count) triples, so every node must exist before they can be resolved, and
    // the arc matrices are sized against the FULL node count rather than the
    // count so far. THE CLASS IS KEPT: `TransitionParam` stores an
    // (nnodes x nclasses) matrix per mode, as MATLAB's `enablingConditions{m}`
    // does, so a mode requiring two Class1 tokens at a place is not satisfied by
    // Class2 tokens sitting there. An arc naming a class the model does not
    // declare is an error rather than a silent drop.
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (node_type[i] != "Transition") continue;
        const json& nd = nodes[i];
        const std::string name = nd.at("name").get<std::string>();
        if (!nd.contains("modes") || nd.at("modes").empty())
            throw InputError("network_reader: transition '" + name + "' declares no mode");
        const json& modes = nd.at("modes");
        qn::TransitionParam<T> tp;
        tp.nmodes = modes.size();
        const std::size_t nn = nodes.size();
        const double inf = std::numeric_limits<double>::infinity();
        for (std::size_t m = 0; m < modes.size(); ++m) {
            const json& mj = modes[m];
            tp.modenames.push_back(mj.value("name", std::string("Mode") + std::to_string(m + 1)));
            // A marking-dependent firing rate g_m(marking) crosses the wire as
            // the MATERIALIZED lattice `{slots, cutoffs, scaling}`, because a
            // handle has no JSON form: `slots` names the enabling (place,class)
            // pairs that g reads, `cutoffs` caps each of them, and `scaling` is
            // keyed by the comma-joined 0-based counts in slot order. Rebuilt
            // here into the same closure MATLAB's `firingdep_table_to_handle`
            // and the Python `_firingdep_table_to_handle` build, over the
            // node-indexed marking `state_events.h` passes in.
            tp.firingdep.push_back(std::function<T(const std::vector<T>&)>());
            if (mj.contains("firingRateDependence")) {
                const json& frm = mj.at("firingRateDependence");
                if (!frm.contains("slots") || !frm.contains("scaling"))
                    throw InputError("network_reader: the marking-dependent firing rate of mode " +
                                     std::to_string(m + 1) + " of transition '" + name +
                                     "' carries no 'slots'/'scaling' lattice");
                std::vector<std::size_t> slot_node;
                for (const json& sm : frm.at("slots")) {
                    const std::string on = sm.at("node").get<std::string>();
                    const std::map<std::string, std::size_t>::const_iterator ni = node_idx.find(on);
                    if (ni == node_idx.end())
                        throw InputError("network_reader: the marking-dependent firing rate of "
                                         "transition '" +
                                         name + "' reads node '" + on +
                                         "', which the model does not declare");
                    slot_node.push_back(ni->second - 1);
                }
                std::vector<long> cutoffs;
                if (frm.contains("cutoffs"))
                    for (const json& c : frm.at("cutoffs")) cutoffs.push_back(c.get<long>());
                std::map<std::string, double> table;
                const json& sc = frm.at("scaling");
                for (json::const_iterator it = sc.begin(); it != sc.end(); ++it)
                    table[it.key()] = it.value().get<double>();
                tp.firingdep.back() = [slot_node, cutoffs,
                                       table](const std::vector<T>& mk) -> T {
                    std::string key;
                    for (std::size_t s = 0; s < slot_node.size(); ++s) {
                        long c = slot_node[s] < mk.size()
                                     ? static_cast<long>(std::llround(
                                           num_traits<T>::to_double(mk[slot_node[s]])))
                                     : 0L;
                        if (c < 0) c = 0;
                        if (s < cutoffs.size() && c > cutoffs[s]) c = cutoffs[s];
                        if (s) key += ',';
                        key += std::to_string(c);
                    }
                    const std::map<std::string, double>::const_iterator hit = table.find(key);
                    // A marking outside the tabulated box is NEUTRAL, not zero:
                    // the same default the three reference readers apply.
                    return hit == table.end() ? num_traits<T>::from_int(1)
                                              : num_traits<T>::from_double(hit->second);
                };
            }
            const bool immediate =
                mj.value("timingStrategy", std::string("TIMED")) == "IMMEDIATE";
            tp.timing.push_back(immediate ? lang::TimingStrategy::IMMEDIATE
                                          : lang::TimingStrategy::TIMED);
            lang::Distrib<T> proc = lang::Distrib<T>::disabled_dist();
            if (mj.contains("distribution")) proc = detail::dist_from_json<T>(mj.at("distribution"));
            tp.firingproc.push_back(proc);
            tp.firingphases.push_back(immediate || proc.disabled ? 0
                                                                 : lang::dist_to_map(proc).order());
            tp.nmodeservers.push_back(detail::num_value(mj, "numServers", 1.0));
            // ONE, not zero: `Transition.addMode` gives a new mode priority 1 in
            // all four codebases (MATLAB Transition.m:88, python nodes.py:2660,
            // JAR Transition.java:114), and the three readers leave that default
            // in place when the key is absent. Defaulting to 0 here made an
            // omitted key mean a DIFFERENT mode than it means everywhere else,
            // and firing priority selects which immediate mode fires, so the
            // marking process itself changed rather than a reported decimal.
            tp.firingprio.push_back(mj.value("firingPriority", 1.0));
            tp.fireweight.push_back(num_traits<T>::from_double(mj.value("firingWeight", 1.0)));
            const std::size_t K = classes.size();
            Matrix<T> enab(nn, K, num_traits<T>::from_int(0));
            Matrix<T> inhib(nn, K, num_traits<T>::from_double(inf));
            Matrix<T> fire(nn, K, num_traits<T>::from_int(0));
            const char* kArcKey[3] = {"enablingConditions", "inhibitingConditions",
                                      "firingOutcomes"};
            for (int which = 0; which < 3; ++which) {
                if (!mj.contains(kArcKey[which])) continue;
                for (const json& arc : mj.at(kArcKey[which])) {
                    const std::string on = arc.at("node").get<std::string>();
                    const std::map<std::string, std::size_t>::const_iterator ni = node_idx.find(on);
                    if (ni == node_idx.end())
                        throw InputError("network_reader: transition '" + name + "' names node '" +
                                         on + "', which the model does not declare");
                    // AN ARC WITHOUT A CLASS IS THE FIRST CLASS'S, which is what
                    // a single-class document means by omitting it; a name the
                    // model never declared is an error, because dropping the arc
                    // would build a net with one fewer precondition.
                    std::size_t rr = 1;
                    if (arc.contains("class")) {
                        const std::string cn = arc.at("class").get<std::string>();
                        const std::map<std::string, std::size_t>::const_iterator ci =
                            class_idx.find(cn);
                        if (ci == class_idx.end())
                            throw InputError("network_reader: transition '" + name +
                                             "' names class '" + cn +
                                             "', which the model does not declare");
                        rr = ci->second;
                    }
                    const T cnt = num_traits<T>::from_double(arc.value("count", 1.0));
                    const std::size_t q = ni->second - 1;
                    if (which == 0) enab(q, rr - 1) = enab(q, rr - 1) + cnt;
                    else if (which == 1) inhib(q, rr - 1) = cnt;  // a threshold, not a count
                    else fire(q, rr - 1) = fire(q, rr - 1) + cnt;
                }
            }
            tp.enabling.push_back(enab);
            tp.inhibiting.push_back(inhib);
            tp.firing.push_back(fire);
        }
        node_idx[name] = net.add_transition(name, tp);
    }

    // -- Pass 2: per-node parameters (service/arrival, servers, CS matrix). ---
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const json& nd = nodes[i];
        const std::string& type = node_type[i];
        const std::size_t idx = node_idx.at(nd.at("name").get<std::string>());

        if (nd.contains("servers") && (type == "Queue"))
            net.set_number_of_servers(idx, detail::num_from_json(nd.at("servers")));

        if (nd.contains("service")) {
            const json& svc = nd.at("service");
            for (auto it = svc.begin(); it != svc.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit == class_idx.end())
                    throw InputError("network_reader: service names unknown class '" + it.key() +
                                     "' at node '" + nd.at("name").get<std::string>() + "'");
                const lang::Distrib<T> d = detail::dist_from_json<T>(it.value());
                if (type == "Source") net.set_arrival(idx, cit->second, d);
                else net.set_service(idx, cit->second, d);
            }
        }

        if (type == "Queue" && nd.contains("scheduling") &&
            nd.at("scheduling").get<std::string>() == "POLLING") {
            if (nd.contains("pollingType"))
                net.set_polling_type(idx, detail::polling_from_json(nd.at("pollingType").get<std::string>()),
                                     nd.value("pollingPar", 0));
        }

        // Per-class scheduling weights (DPS / GPS). Without these the AMVA DPS
        // path either errors (no weights) or degenerates to plain PS.
        if (nd.contains("schedParams") && (type == "Queue" || type == "Delay")) {
            const json& sp = nd.at("schedParams");
            for (auto it = sp.begin(); it != sp.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_sched_param(idx, cit->second,
                                        num_traits<T>::from_double(it.value().get<double>()));
            }
        }

        // Finite-buffer blocking. `dropRule` per class (blockingAfterService ->
        // BAS makes mva_is_bas_model true, routing the model to solver_sqd) and
        // `buffer` the total station capacity.
        // The writers emit these for `isa(node,'Station')`, which is wider than
        // Queue and Delay: a Join and a queueing Place are stations too, and a
        // rule read only at a Queue is a rule silently dropped at those.
        const bool station_node =
            type == "Queue" || type == "Delay" || type == "Join" || type == "Place";
        if (nd.contains("dropRule") && station_node) {
            const json& dr = nd.at("dropRule");
            for (auto it = dr.begin(); it != dr.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_drop_rule(idx, cit->second,
                                      detail::drop_from_json(it.value().get<std::string>()));
            }
        }
        if (nd.contains("buffer") && station_node)
            net.set_capacity(idx, detail::num_from_json(nd.at("buffer")));
        // Per-class buffer capacity, the refinement of `buffer`.
        if (nd.contains("classCap")) {
            const json& cc = nd.at("classCap");
            for (auto it = cc.begin(); it != cc.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_class_capacity(idx, cit->second, detail::num_from_json(it.value()));
            }
        }

        // -- Impatience, balking and the retrial orbit -----------------------
        //
        // Four DIFFERENT populations, and the reference keeps them apart: a
        // patience timer abandons a job that has JOINED the queue, balking
        // refuses to join at all, an orbit impatience abandons a retrial orbit
        // (never a buffer slot), and a batch rejection drops a whole arriving
        // batch. Folding any two together changes which jobs are counted.
        if (nd.contains("patience")) {
            const json& pt = nd.at("patience");
            for (auto it = pt.begin(); it != pt.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit == class_idx.end()) continue;
                const json& pj = it.value();
                const lang::ImpatienceType kind =
                    pj.contains("impatienceType")
                        ? detail::impatience_from_json(pj.at("impatienceType").get<std::string>())
                        : lang::ImpatienceType::RENEGING;
                net.set_patience(idx, cit->second, detail::dist_from_json<T>(pj.at("distribution")),
                                 kind);
            }
        }
        if (nd.contains("orbitImpatience")) {
            const json& oi = nd.at("orbitImpatience");
            for (auto it = oi.begin(); it != oi.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_orbit_impatience(idx, cit->second, detail::dist_from_json<T>(it.value()));
            }
        }
        if (nd.contains("batchRejectProb")) {
            const json& br = nd.at("batchRejectProb");
            for (auto it = br.begin(); it != br.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_batch_reject(idx, cit->second,
                                         num_traits<T>::from_double(it.value().get<double>()));
            }
        }
        if (nd.contains("balking")) {
            const json& bk = nd.at("balking");
            for (auto it = bk.begin(); it != bk.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit == class_idx.end()) continue;
                const json& bj = it.value();
                std::vector<typename qn::Station<T>::BalkingThreshold> ths;
                if (bj.contains("thresholds"))
                    for (const json& tj : bj.at("thresholds")) {
                        typename qn::Station<T>::BalkingThreshold th;
                        th.min_jobs = tj.value("minJobs", 0.0);
                        // -1 IS THE WIRE'S INFINITY, not a count: the writer maps
                        // an unbounded upper end to it rather than emitting Inf,
                        // which JSON has no literal for.
                        th.max_jobs = tj.value("maxJobs", -1.0);
                        th.probability = num_traits<T>::from_double(tj.value("probability", 1.0));
                        ths.push_back(th);
                    }
                net.set_balking(idx, cit->second,
                                detail::balking_from_json(bj.at("strategy").get<std::string>()), ths);
            }
        }
        if (nd.contains("retrial")) {
            const json& rt = nd.at("retrial");
            for (auto it = rt.begin(); it != rt.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit == class_idx.end()) continue;
                const lang::Distrib<T> delay = detail::dist_from_json<T>(it.value().at("delay"));
                // `sn.retrialMu` is the RATE of that delay, which is what the
                // state layer multiplies the orbit population by.
                const T rate = T(num_traits<T>::from_int(1) / delay.mean);
                net.set_retrial(idx, cit->second, delay, rate,
                                int(detail::num_value(it.value(), "maxAttempts", 0.0)));
            }
        }

        // -- Setup / delay-off, switchover and heterogeneous servers ---------
        if (nd.contains("setupTime")) {
            if (!nd.contains("delayOffTime"))
                throw InputError(
                    "network_reader: node '" + nd.at("name").get<std::string>() +
                    "' declares a setup time with no delay-off time; a server that never powers "
                    "down never pays the setup, so the pair is meaningless alone");
            const json& su = nd.at("setupTime");
            const json& doff = nd.at("delayOffTime");
            for (auto it = su.begin(); it != su.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit == class_idx.end() || !doff.contains(it.key())) continue;
                net.set_setup_delayoff(idx, cit->second, detail::dist_from_json<T>(it.value()),
                                       detail::dist_from_json<T>(doff.at(it.key())));
            }
        }
        // -- Server breakdown / repair ---------------------------------------
        // The wire form is `linemodel_save`'s: {failure, repair, downService?},
        // with `downService` keyed by CLASS NAME. MATLAB writes it per class
        // even when `setBreakdown` was given one distribution for every class,
        // flattening the class-independent form, so there is only one shape to
        // read here.
        if (nd.contains("breakdown")) {
            const json& bd = nd.at("breakdown");
            if (!bd.contains("failure") || !bd.contains("repair"))
                throw InputError(
                    "network_reader: node '" + nd.at("name").get<std::string>() +
                    "' declares a breakdown without both a failure and a repair time; a server "
                    "that never recovers is an absorbing model, not a breakdown");
            std::vector<lang::Distrib<T> > down(class_idx.size(),
                                                lang::Distrib<T>::disabled_dist());
            if (bd.contains("downService")) {
                const json& ds = bd.at("downService");
                for (auto it = ds.begin(); it != ds.end(); ++it) {
                    auto cit = class_idx.find(it.key());
                    if (cit == class_idx.end()) continue;
                    if (cit->second >= 1 && cit->second <= down.size())
                        down[cit->second - 1] = detail::dist_from_json<T>(it.value());
                }
            }
            net.set_breakdown(idx, detail::dist_from_json<T>(bd.at("failure")),
                              detail::dist_from_json<T>(bd.at("repair")), down);
        }
        if (nd.contains("switchoverTimes")) {
            // A SWITCHOVER IS A POLLING SERVER'S WALK, AND NOTHING ELSE READS
            // ONE. Every consumer in this port -- the state machinery, the SSA
            // events, the JMT writer -- reaches it through the polling buffers,
            // and so does the reference: `State.afterEventStation` handles
            // SWITCH only inside `pollingInfo`, and `writeJSIM` warns and drops
            // the times on an ordinary Server. Declared away from a POLLING
            // queue the block is therefore inert in both codebases, and it is
            // dropped here with the reference's own warning rather than
            // refused: refusing made switchover_basic, which the reference
            // solves and reports, unsolvable in C++.
            const bool polling = type == "Queue" &&
                                 nd.value("scheduling", std::string()) == "POLLING";
            if (!polling) {
                std::cerr << "[LINE] Warning: node '" << nd.at("name").get<std::string>()
                          << "' declares switchover times but is not POLLING-scheduled; a "
                          << "switchover is the walk between a polling server's buffers, so "
                          << "the times are ignored." << std::endl;
            } else {
                // Indexed by the class the server LEAVES. The non-polling form
                // also names the class it moves TO, which this port's per-class
                // vector cannot express, so a genuinely (from, to)-dependent
                // walk is refused rather than collapsed onto its `from` alone.
                std::map<std::size_t, lang::Distrib<T> > walk;
                std::map<std::size_t, std::string> first_to;
                for (const json& so : nd.at("switchoverTimes")) {
                    auto cit = class_idx.find(so.at("from").get<std::string>());
                    if (cit == class_idx.end()) continue;
                    const std::string to = so.value("to", std::string());
                    if (walk.count(cit->second) && first_to[cit->second] != to)
                        throw UnsupportedError(
                            "network_reader: node '" + nd.at("name").get<std::string>() +
                            "' declares a switchover that depends on the class moved TO as well "
                            "as the one moved FROM; this port carries one walk per departing "
                            "class");
                    walk[cit->second] = detail::dist_from_json<T>(so.at("distribution"));
                    first_to[cit->second] = to;
                }
                for (const auto& kv : walk) net.set_switchover(idx, kv.first, kv.second);
            }
        }
        if (nd.contains("serverTypes")) {
            for (const json& st : nd.at("serverTypes")) {
                typename qn::Station<T>::ServerType stype;
                stype.name = st.value("name", std::string());
                stype.count = st.value("count", 1.0);
                if (st.contains("compatibleClasses")) {
                    stype.compatible.assign(classes.size(), false);
                    for (const json& cn : st.at("compatibleClasses")) {
                        auto cit = class_idx.find(cn.get<std::string>());
                        if (cit != class_idx.end()) stype.compatible[cit->second - 1] = true;
                    }
                }
                if (st.contains("service")) {
                    stype.service.assign(classes.size(), lang::Distrib<T>::disabled_dist());
                    const json& sv = st.at("service");
                    for (auto it = sv.begin(); it != sv.end(); ++it) {
                        auto cit = class_idx.find(it.key());
                        if (cit != class_idx.end())
                            stype.service[cit->second - 1] = detail::dist_from_json<T>(it.value());
                    }
                }
                net.add_server_type(idx, stype);
            }
            if (nd.contains("heteroSchedPolicy"))
                net.set_hetero_sched_policy(
                    idx, detail::hetero_from_json(nd.at("heteroSchedPolicy").get<std::string>()));
        }
        // Job parallelism: servers seized at once by a job, per class. Read after
        // the pools, which size the station and so bound the admissible values.
        if (nd.contains("serverParallelism")) {
            const json& sp = nd.at("serverParallelism");
            for (auto it = sp.begin(); it != sp.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_server_parallelism(idx, cit->second, it.value().get<std::size_t>());
            }
        }

        // -- The Source refinements: batch size and the MMAP mark binding ----
        if (nd.contains("arrivalBatch")) {
            const json& ab = nd.at("arrivalBatch");
            for (auto it = ab.begin(); it != ab.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_arrival_batch(idx, cit->second, detail::dist_from_json<T>(it.value()));
            }
        }
        if (nd.contains("markedClasses")) {
            std::vector<std::size_t> marks;
            for (const json& cn : nd.at("markedClasses")) {
                auto cit = class_idx.find(cn.get<std::string>());
                if (cit == class_idx.end())
                    throw InputError("network_reader: node '" + nd.at("name").get<std::string>() +
                                     "' binds mark " + std::to_string(marks.size() + 1) +
                                     " to class '" + cn.get<std::string>() +
                                     "', which the model does not declare");
                marks.push_back(cit->second);
            }
            net.set_marked_classes(idx, marks);
        }

        // -- The order-independent / pass-and-swap station -------------------
        //
        // mu(c) crosses as a macrostate TABLE, so `oiCutoffs` is not decoration:
        // it is the lattice the table was materialized over and the point beyond
        // which the rate saturates.
        if (nd.contains("oiServiceRate")) {
            std::vector<int> cut(classes.size(), 10);
            if (nd.contains("oiCutoffs")) {
                const std::vector<double> cv = detail::num_vec_from_json<double>(nd.at("oiCutoffs"));
                for (std::size_t r = 0; r < cut.size() && r < cv.size(); ++r)
                    cut[r] = static_cast<int>(std::lround(cv[r]));
            }
            std::vector<std::vector<bool> > swap;
            if (nd.contains("swapGraph")) {
                const Matrix<T> sg = detail::mat_from_json<T>(nd.at("swapGraph"));
                swap.assign(sg.rows(), std::vector<bool>(sg.cols(), false));
                for (std::size_t a = 0; a < sg.rows(); ++a)
                    for (std::size_t b = 0; b < sg.cols(); ++b)
                        swap[a][b] = num_traits<T>::to_double(sg(a, b)) != 0.0;
            }
            net.set_pas(idx,
                        detail::oi_rate_from_json<T>(nd.at("oiServiceRate"), cut, classes.size()),
                        swap);
        }

        // -- The declared state -----------------------------------------------
        if (nd.contains("initialState") && type == "Place") {
            net.set_initial_marking(idx, detail::num_vec_from_json<T>(nd.at("initialState")));
        } else if (nd.contains("initialState") && !nd.contains("stateSpace")) {
            // A stateful node that is not a Place carries its declared state in
            // the same key, and this struct spells a declared state as a
            // ONE-ROW state space under a prior of one (`sn_state.h`, the CTMC
            // transient's own reading). MATLAB, the JAR and Python all write
            // the key for every stateful node, so ignoring it outside a Place
            // dropped the initialization of every station: a model saved after
            // initFromMarginal came back here as the default one, with the
            // difference visible only in the numbers.
            const std::vector<T> row = detail::num_vec_from_json<T>(nd.at("initialState"));
            if (!row.empty()) {
                Matrix<T> space(1, row.size());
                for (std::size_t k = 0; k < row.size(); ++k) space(0, k) = row[k];
                net.set_state_prior(idx, space, std::vector<T>(1, num_traits<T>::from_double(1.0)));
            }
        }
        if (nd.contains("statePrior") || nd.contains("stateSpace")) {
            if (!nd.contains("statePrior") || !nd.contains("stateSpace"))
                throw InputError(
                    "network_reader: node '" + nd.at("name").get<std::string>() +
                    "' declares one of stateSpace / statePrior without the other; the prior is a "
                    "distribution over the ROWS of that space and means nothing alone");
            net.set_state_prior(idx, detail::mat_from_json<T>(nd.at("stateSpace")),
                                detail::num_vec_from_json<T>(nd.at("statePrior")));
        }
        if (nd.contains("departureDiscipline")) {
            const json& dd = nd.at("departureDiscipline");
            for (auto it = dd.begin(); it != dd.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end())
                    net.set_departure_discipline(
                        idx, cit->second, detail::departure_from_json(it.value().get<std::string>()));
            }
        }

        // Rate scaling. `loadDependence` carries the lldscaling vector itself,
        // indexed from population one; `classDependence` and `jointDependence`
        // carry beta_r(n) and eta_i(n) materialized over the per-class box
        // lattice, because a function handle cannot cross JSON.
        for (int kind = 0; kind < 3; ++kind) {
            static const char* kKeys[] = {"loadDependence", "classDependence",
                                          "jointDependence"};
            static const char* kTypes[] = {"loadDependent", "classDependent",
                                           "jointDependent"};
            if (!nd.contains(kKeys[kind])) continue;
            if (type != "Queue" && type != "Delay")
                throw UnsupportedError(
                    std::string("network_reader: node '") + nd.at("name").get<std::string>() +
                    "' carries '" + kKeys[kind] +
                    "', which scales a SERVICE rate and is only meaningful at a Queue or a "
                    "Delay");
            const detail::json& blk = nd.at(kKeys[kind]);
            if (blk.value("type", std::string(kTypes[kind])) != kTypes[kind])
                throw UnsupportedError(std::string("network_reader: '") + kKeys[kind] +
                                       "' at node '" + nd.at("name").get<std::string>() +
                                       "' declares type '" +
                                       blk.value("type", std::string("?")) +
                                       "', which this reader does not implement");
            if (!blk.contains("scaling"))
                throw InputError(std::string("network_reader: '") + kKeys[kind] +
                                 "' at node '" + nd.at("name").get<std::string>() +
                                 "' carries no 'scaling'");
            if (kind == 0) {
                net.set_load_dependence(idx, detail::num_vec_from_json<T>(blk.at("scaling")));
                continue;
            }
            // An absent `cutoffs` leaves the lattice at one job per class, which
            // is what the Python reader assumes for the same legacy JSON.
            std::vector<int> cut(classes.size(), 1);
            if (blk.contains("cutoffs")) {
                const std::vector<double> cv = detail::num_vec_from_json<double>(blk.at("cutoffs"));
                for (std::size_t r = 0; r < cut.size() && r < cv.size(); ++r)
                    cut[r] = static_cast<int>(std::lround(cv[r]));
            }
            const lang::CdScaling<T> fun =
                detail::cd_scaling_from_json<T>(blk.at("scaling"), cut, classes.size());
            const std::vector<T> peak = detail::cd_peak_from_json<T>(blk, blk.at("scaling"));
            if (kind == 1) net.set_class_dependence(idx, fun, peak);
            else net.set_joint_dependence(idx, fun, peak);
        }
        // Immediate feedback, a per-class map on the node. The writers emit it
        // for a Queue only, and only for the classes it is set on, so an absent
        // entry means "not set" rather than false.
        if (nd.contains("immediateFeedback")) {
            const json& imf = nd.at("immediateFeedback");
            for (auto it = imf.begin(); it != imf.end(); ++it) {
                auto cit = class_idx.find(it.key());
                if (cit != class_idx.end() && it.value().get<bool>())
                    net.set_immediate_feedback(idx, cit->second);
            }
        }
    }

    // -- Pass 3: routing. ----------------------------------------------------
    const json& routing = model.at("routing");
    qn::RoutingMatrix<T> P;
    const json& mat = routing.at("matrix");
    for (auto ck = mat.begin(); ck != mat.end(); ++ck) {
        // Key "SrcClass,DstClass"; a bare "Class" means same-class routing.
        const std::string key = ck.key();
        const std::size_t comma = key.find(',');
        const std::string cs = comma == std::string::npos ? key : key.substr(0, comma);
        const std::string cd = comma == std::string::npos ? key : key.substr(comma + 1);
        const std::size_t r = class_idx.at(cs);
        const std::size_t s = class_idx.at(cd);
        for (auto si = ck.value().begin(); si != ck.value().end(); ++si) {
            const std::size_t from = node_idx.at(si.key());
            for (auto di = si.value().begin(); di != si.value().end(); ++di)
                P.set(r, s, from, node_idx.at(di.key()),
                      num_traits<T>::from_double(di.value().get<double>()));
        }
    }
    net.link(P);

    // -- Pass 3b: the non-probabilistic dispatchers. -------------------------
    //
    // `routing` above carries the MATRIX, which is what a PROB dispatcher is
    // fully described by. A RROBIN, JSQ or power-of-d dispatcher has the same
    // matrix and a different rule, so reading only the matrix turns every one
    // of them into probabilistic routing with no diagnostic. `link()` has
    // already run, so these overwrite the per-class strategy it defaulted to.
    if (model.contains("routingStrategies")) {
        const json& rs = model.at("routingStrategies");
        for (auto ni = rs.begin(); ni != rs.end(); ++ni) {
            auto nit = node_idx.find(ni.key());
            for (auto ci = ni.value().begin(); ci != ni.value().end(); ++ci) {
                const lang::RoutingStrategy rst =
                    detail::routing_from_json(ci.value().get<std::string>());
                // DISABLED IS DERIVED, not declared: it marks a (node, class)
                // pair the class never visits, the refresh re-derives it, and
                // every reference loader skips it. Writers before BUG-94 also
                // emitted it for the AUTO-ADDED class-switch nodes they do not
                // put in `nodes`, so an entry naming an unknown node is only
                // tolerated for this one value -- anything else is still a
                // dangling reference and is named.
                if (nit == node_idx.end()) {
                    if (rst == lang::RoutingStrategy::DISABLED) continue;
                    throw InputError("network_reader: routingStrategies names node '" + ni.key() +
                                     "', which the model does not declare");
                }
                if (rst == lang::RoutingStrategy::DISABLED) continue;
                auto cit = class_idx.find(ci.key());
                if (cit == class_idx.end()) continue;
                net.set_routing(nit->second, cit->second, rst);
            }
        }
    }
    if (model.contains("routingWeights")) {
        const json& rw = model.at("routingWeights");
        for (auto ni = rw.begin(); ni != rw.end(); ++ni) {
            auto nit = node_idx.find(ni.key());
            if (nit == node_idx.end()) continue;
            for (auto ci = ni.value().begin(); ci != ni.value().end(); ++ci) {
                auto cit = class_idx.find(ci.key());
                if (cit == class_idx.end()) continue;
                std::map<std::size_t, double> w;
                for (auto di = ci.value().begin(); di != ci.value().end(); ++di) {
                    auto dit = node_idx.find(di.key());
                    if (dit != node_idx.end()) w[dit->second] = di.value().get<double>();
                }
                net.set_routing_weights(nit->second, cit->second, w);
            }
        }
    }
    if (model.contains("routingParams")) {
        const json& rp = model.at("routingParams");
        for (auto ni = rp.begin(); ni != rp.end(); ++ni) {
            auto nit = node_idx.find(ni.key());
            if (nit == node_idx.end()) continue;
            for (auto ci = ni.value().begin(); ci != ni.value().end(); ++ci) {
                auto cit = class_idx.find(ci.key());
                if (cit == class_idx.end() || !ci.value().contains("d")) continue;
                net.set_routing_param(nit->second, cit->second, ci.value().at("d").get<int>());
            }
        }
    }

    // -- Pass 3b3: Krzesinski state-dependent routing. -----------------------
    //
    // Restored AFTER `link()` and after the dispatchers above: `link` writes the
    // uniform placeholder into the entry row, which the declaration supersedes,
    // and `routingStrategies` has already named that row SDR without saying what
    // the subnetwork looks like. Every centre travels by NODE NAME, so a node
    // reordering on the writing side cannot shift one, and branch index 1 is the
    // complement M-V and arrives as an empty list.
    // See _kb/16-state-dependent-routing.md
    if (model.contains("stateDepRouting")) {
        const json& sd = model.at("stateDepRouting");
        for (const char* k : {"entry", "departure", "class", "branches", "level", "C", "d"})
            if (!sd.contains(k))
                throw InputError(std::string("network_reader: 'stateDepRouting' carries no '") + k +
                                 "'");
        auto node_of = [&](const std::string& nm) -> std::size_t {
            const std::map<std::string, std::size_t>::const_iterator it = node_idx.find(nm);
            if (it == node_idx.end())
                throw InputError("network_reader: 'stateDepRouting' names node '" + nm +
                                 "', which the model does not declare");
            return it->second;
        };
        const std::map<std::string, std::size_t>::const_iterator ci =
            class_idx.find(sd.at("class").get<std::string>());
        if (ci == class_idx.end())
            throw InputError("network_reader: 'stateDepRouting' names class '" +
                             sd.at("class").get<std::string>() +
                             "', which the model does not declare");
        std::vector<std::vector<std::size_t>> branches;
        for (const json& b : sd.at("branches")) {
            std::vector<std::size_t> centres;
            for (const json& nm : b) centres.push_back(node_of(nm.get<std::string>()));
            branches.push_back(centres);
        }
        std::vector<std::size_t> level;
        for (const json& v : sd.at("level")) level.push_back(static_cast<std::size_t>(v.get<double>()));
        std::vector<double> C;
        for (const json& v : sd.at("C")) C.push_back(v.get<double>());
        const json& dj = sd.at("d");
        Matrix<double> d(dj.size(), dj.empty() ? 0 : dj.at(0).size(), 0.0);
        for (std::size_t t = 0; t < dj.size(); ++t) {
            if (dj.at(t).size() != d.cols())
                throw InputError("network_reader: 'stateDepRouting' carries a ragged coefficient "
                                 "matrix d");
            for (std::size_t b = 0; b < d.cols(); ++b) d(t, b) = dj.at(t).at(b).get<double>();
        }
        net.set_state_dep_routing(node_of(sd.at("entry").get<std::string>()),
                                  node_of(sd.at("departure").get<std::string>()), branches, level, C,
                                  d, ci->second);
    }

    // -- Pass 3b2: the global (Whittle) dependence. --------------------------
    //
    // phi(n) is rebuilt from the materialized slot lattice written by
    // `network_to_json`. The slots carry station and class NAMES, resolved
    // through this model's own index spaces, so a reordering on the writing side
    // cannot silently shift a coordinate. The population is clamped to the
    // tabulated cutoffs, which is the same saturation the writer's box lattice
    // declares.
    if (model.contains("globalDependence")) {
        const json& blk = model.at("globalDependence");
        if (blk.value("type", std::string("globalDependent")) != "globalDependent")
            throw UnsupportedError(
                std::string("network_reader: 'globalDependence' declares type '") +
                blk.value("type", std::string("?")) + "', which this reader does not implement");
        if (!blk.contains("scaling"))
            throw InputError("network_reader: 'globalDependence' carries no 'scaling'");
        const qn::NetworkStruct<T>& sn0 = net.get_struct();
        const std::size_t M = sn0.nstations, K = sn0.nclasses;
        std::vector<std::size_t> slot_st, slot_cl;
        if (blk.contains("slots"))
            for (const json& sm : blk.at("slots")) {
                auto nit = node_idx.find(sm.at("station").get<std::string>());
                auto cit = class_idx.find(sm.at("class").get<std::string>());
                if (nit == node_idx.end() || cit == class_idx.end())
                    throw InputError(
                        "network_reader: 'globalDependence' names a station or class the model "
                        "does not declare");
                slot_st.push_back(sn0.nodes[nit->second - 1].station - 1);
                slot_cl.push_back(cit->second - 1);
            }
        const std::size_t P = slot_st.size();
        std::vector<int> cuts(P, 0);
        if (blk.contains("cutoffs")) {
            const std::vector<double> cv = detail::num_vec_from_json<double>(blk.at("cutoffs"));
            for (std::size_t d = 0; d < P && d < cv.size(); ++d)
                cuts[d] = static_cast<int>(std::lround(cv[d]));
        }
        const int wcut = blk.value("cutoff", 10);
        std::map<std::string, std::vector<T>> tbl;
        for (auto it = blk.at("scaling").begin(); it != blk.at("scaling").end(); ++it)
            tbl[it.key()] = detail::num_vec_from_json<T>(it.value());
        std::vector<T> peak(M * K, num_traits<T>::from_int(1));
        if (blk.contains("peak")) {
            const std::vector<T> pv = detail::num_vec_from_json<T>(blk.at("peak"));
            for (std::size_t j = 0; j < peak.size() && j < pv.size(); ++j) peak[j] = pv[j];
        }
        const std::vector<T> ones(M * K, num_traits<T>::from_int(1));
        net.set_global_dependence(
            [slot_st, slot_cl, cuts, tbl, ones, K, P](const std::vector<T>& n) {
                std::string key;
                if (P == 0) key = "0";
                else
                    for (std::size_t d = 0; d < P; ++d) {
                        long x = std::lround(num_traits<T>::to_double(n[slot_st[d] * K + slot_cl[d]]));
                        if (x < 0) x = 0;
                        if (x > cuts[d]) x = cuts[d];
                        if (d) key += ',';
                        key += std::to_string(x);
                    }
                typename std::map<std::string, std::vector<T>>::const_iterator it = tbl.find(key);
                return it == tbl.end() ? ones : it->second;
            },
            peak, wcut);
    }

    // -- Pass 3c: the finite capacity regions. -------------------------------
    //
    // A region caps the jobs held ACROSS a set of stations, which no per-station
    // capacity expresses. Dropping the block is the failure this reader's key
    // gate was written for: an FCR model then solves as the UNBOUNDED one and
    // reports a confident, wrong queue length.
    if (model.contains("finiteCapacityRegions")) {
        for (const json& rj : model.at("finiteCapacityRegions")) {
            std::vector<std::size_t> members;
            // The per-station `classCap` refines the region's own cap at that
            // station; the struct carries ONE cap per (station, class), so the
            // tightest of the two is what the region actually enforces.
            std::vector<double> cap(classes.size(), -1.0);
            auto read_class_map = [&](const json& blk, std::vector<double>& dst) {
                for (auto it = blk.begin(); it != blk.end(); ++it) {
                    auto cit = class_idx.find(it.key());
                    if (cit == class_idx.end()) continue;
                    const double v = it.value().get<double>();
                    double& slot = dst[cit->second - 1];
                    slot = slot == -1.0 ? v : std::min(slot, v);
                }
            };
            if (rj.contains("classMaxJobs")) read_class_map(rj.at("classMaxJobs"), cap);
            for (const json& sj : rj.at("stations")) {
                auto nit = node_idx.find(sj.at("node").get<std::string>());
                if (nit == node_idx.end())
                    throw InputError("network_reader: finite capacity region '" +
                                     rj.value("name", std::string("?")) + "' names node '" +
                                     sj.at("node").get<std::string>() +
                                     "', which the model does not declare");
                members.push_back(nit->second);
                if (sj.contains("classCap")) read_class_map(sj.at("classCap"), cap);
            }
            std::vector<double> mem(classes.size(), -1.0);
            if (rj.contains("classMaxMemory")) read_class_map(rj.at("classMaxMemory"), mem);
            std::vector<T> size(classes.size(), num_traits<T>::from_int(1));
            std::vector<T> weight(classes.size(), num_traits<T>::from_int(1));
            for (const json& sj : rj.at("stations")) {
                if (sj.contains("classSize"))
                    for (auto it = sj.at("classSize").begin(); it != sj.at("classSize").end(); ++it) {
                        auto cit = class_idx.find(it.key());
                        if (cit != class_idx.end())
                            size[cit->second - 1] =
                                num_traits<T>::from_double(it.value().get<double>());
                    }
                if (sj.contains("classWeight"))
                    for (auto it = sj.at("classWeight").begin(); it != sj.at("classWeight").end();
                         ++it) {
                        auto cit = class_idx.find(it.key());
                        if (cit != class_idx.end())
                            weight[cit->second - 1] =
                                num_traits<T>::from_double(it.value().get<double>());
                    }
            }
            std::vector<lang::DropStrategy> rule(classes.size(), lang::DropStrategy::WAITQ);
            if (rj.contains("dropRule"))
                for (auto it = rj.at("dropRule").begin(); it != rj.at("dropRule").end(); ++it) {
                    auto cit = class_idx.find(it.key());
                    if (cit != class_idx.end())
                        rule[cit->second - 1] = detail::drop_from_json(it.value().get<std::string>());
                }
            const std::size_t reg =
                net.add_region(members, cap, rj.value("globalMaxJobs", -1.0), rule, mem, size,
                               rj.value("globalMaxMemory", -1.0),
                               rj.value("name", std::string()));
            net.set_region_weights(reg, weight);
            if (rj.contains("constraintA") && rj.contains("constraintB"))
                net.set_region_constraint(reg, detail::mat_from_json<T>(rj.at("constraintA")),
                                          detail::num_vec_from_json<T>(rj.at("constraintB")));
        }
    }

    // -- Pass 4: the declared rewards. ---------------------------------------
    //
    // THE DECLARATIVE FORM ONLY, `{name, type, node, class?}`, which is all the
    // writers emit: a reward defined from a bare lambda has no reproducible form
    // and `linemodel_save` warns and omits it rather than writing something that
    // would be wrong on reload. Each template is rebuilt here as a function of
    // the AGGREGATE state row -- the per-(station, class) counts in
    // `(ist-1)*K + k` order that `set_reward` is defined over -- so the value is
    // computed from the chain's own states and not from a mean.
    if (model.contains("rewards")) {
        const std::size_t K = classes.size();
        for (const json& rw : model.at("rewards")) {
            const std::string nm = rw.at("name").get<std::string>();
            const std::string kind = rw.at("type").get<std::string>();
            const std::string node_name = rw.at("node").get<std::string>();
            auto nit = node_idx.find(node_name);
            if (nit == node_idx.end())
                throw InputError("network_reader: reward '" + nm + "' names node '" + node_name +
                                 "', which the model does not declare");
            // The struct carries station -> node and no inverse, so the station
            // is found by scanning it; a node with no station is refused below.
            std::size_t ist = 0;
            for (std::size_t s = 0; s < net.get_struct().station_to_node.size(); ++s)
                if (net.get_struct().station_to_node[s] == nit->second) ist = s + 1;
            if (ist == 0)
                throw UnsupportedError("network_reader: reward '" + nm + "' is declared at node '" +
                                       node_name +
                                       "', which is not a station and so has no job count");
            // A reward with no class covers EVERY class at the station; with one,
            // exactly that class. The two are different quantities, so the
            // absence of the key is carried through rather than defaulted.
            std::size_t cls = 0;  // 0 = all classes
            if (rw.contains("class")) {
                auto cit = class_idx.find(rw.at("class").get<std::string>());
                if (cit == class_idx.end())
                    throw InputError("network_reader: reward '" + nm + "' names class '" +
                                     rw.at("class").get<std::string>() +
                                     "', which the model does not declare");
                cls = cit->second;
            }
            const double nservers = net.get_struct().stations[ist - 1].nservers;
            const double cap = net.get_struct().stations[ist - 1].cap;
            const std::size_t base = (ist - 1) * K;
            // The DECLARATIVE descriptor travels with the lambda so the writer
            // can emit the reward back out; a reward built from a bare function
            // has none and is omitted at save time, as the reference does.
            if (kind == "QLen") {
                net.set_reward(nm, [base, K, cls](const std::vector<T>& n) {
                    if (cls) return n[base + cls - 1];
                    T s = num_traits<T>::from_int(0);
                    for (std::size_t k = 0; k < K; ++k) s = T(s + n[base + k]);
                    return s;
                }, kind, nit->second, cls);
            } else if (kind == "Util") {
                // min(jobs, nservers), the reference's own definition: the number
                // of BUSY servers, not a fraction.
                net.set_reward(nm, [base, K, cls, nservers](const std::vector<T>& n) {
                    T s = num_traits<T>::from_int(0);
                    if (cls) s = n[base + cls - 1];
                    else
                        for (std::size_t k = 0; k < K; ++k) s = T(s + n[base + k]);
                    const T c = num_traits<T>::from_double(nservers);
                    return num_traits<T>::to_double(s) > nservers ? c : s;
                }, kind, nit->second, cls);
            } else if (kind == "Blocking") {
                net.set_reward(nm, [base, K, cap](const std::vector<T>& n) {
                    T s = num_traits<T>::from_int(0);
                    for (std::size_t k = 0; k < K; ++k) s = T(s + n[base + k]);
                    return num_traits<T>::to_double(s) >= cap ? num_traits<T>::from_int(1)
                                                              : num_traits<T>::from_int(0);
                }, kind, nit->second, cls);
            } else {
                throw UnsupportedError(
                    "network_reader: reward '" + nm + "' is of type '" + kind +
                    "', and the reproducible templates are QLen, Util and Blocking; a Custom "
                    "reward wraps an arbitrary function and is not on the wire at all");
            }
        }
    }
    return net;
}

/** Parse a model.json file into a `qn::Network<T>`. */
template <class T>
qn::Network<T> read_network_json(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("network_reader: cannot open " + path);
    detail::json root;
    try {
        in >> root;
    } catch (const detail::json::parse_error& e) {
        throw InputError("network_reader: malformed JSON in " + path + ": " + e.what());
    }
    return build_network_from_json<T>(root);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_NETWORK_READER_H
