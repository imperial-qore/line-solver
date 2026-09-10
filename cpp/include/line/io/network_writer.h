/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_NETWORK_WRITER_H
#define LINE_IO_NETWORK_WRITER_H

/**
 * `qn::NetworkStruct` -> model.json, the inverse of network_reader.h.
 *
 * The wire format is the one `linemodel_save.m`, the Python `save_model` and
 * the JAR `LineModelIO` write, so a model saved here loads in any of the four
 * codebases and a model saved by any of them and read back here is unchanged.
 *
 * WHAT IS AND IS NOT ROUND-TRIPPABLE. Two things on a model are FUNCTIONS, and
 * a function cannot cross JSON:
 *
 *   the rate-scaling handles     beta_r(n), eta_i(n) and the OI mu(c)
 *   a reward built from a lambda
 *
 * The reference solves the first by MATERIALIZING the handle over the per-class
 * box lattice, and this writer does the same, with the same cutoffs (a closed
 * class's population, or 10 for an open one) and the same comma-joined key, so
 * the table is byte-comparable with MATLAB's. The second it solves by writing
 * only the DECLARATIVE rewards and warning about the rest; here a reward with
 * no `kind` is omitted for the same reason -- emitting a guess would produce a
 * model.json that reloads into a different reward.
 *
 * EVERY OTHER KEY IS EXACT. The writer emits exactly the keys the reader
 * consumes, which is what makes the round trip a test rather than a hope: a
 * field added to one side without the other shows up as a read-back mismatch
 * in `test_network_roundtrip.cpp`.
 */

#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/io/network_reader.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace detail {

/** The model.json `type` of a node kind. */
inline const char* node_type_to_json(lang::NodeType t) {
    switch (t) {
        case lang::NodeType::Queue: return "Queue";
        case lang::NodeType::Source: return "Source";
        case lang::NodeType::Delay: return "Delay";
        case lang::NodeType::ClassSwitch: return "ClassSwitch";
        case lang::NodeType::Logger: return "Logger";
        case lang::NodeType::Cache: return "Cache";
        case lang::NodeType::Router: return "Router";
        case lang::NodeType::Fork: return "Fork";
        case lang::NodeType::Place: return "Place";
        case lang::NodeType::Transition: return "Transition";
        case lang::NodeType::Region: return "Region";
        case lang::NodeType::Join: return "Join";
        case lang::NodeType::Sink: return "Sink";
    }
    throw UnsupportedError("network_writer: unnamed node type");
}

inline const char* drop_to_json(lang::DropStrategy d) {
    switch (d) {
        case lang::DropStrategy::WAITQ: return "waitingQueue";
        case lang::DropStrategy::DROP: return "drop";
        case lang::DropStrategy::BAS: return "blockingAfterService";
        case lang::DropStrategy::BBS: return "blockingBeforeService";
        case lang::DropStrategy::RSRD: return "resamplingRepetitiveService";
        case lang::DropStrategy::RETRIAL: return "retrial";
        case lang::DropStrategy::RETRIAL_WITH_LIMIT: return "retrialWithLimit";
    }
    throw UnsupportedError("network_writer: unnamed drop strategy");
}

inline const char* replacement_to_json(lang::ReplacementStrategy r) {
    switch (r) {
        case lang::ReplacementStrategy::RR: return "RR";
        case lang::ReplacementStrategy::FIFO: return "FIFO";
        case lang::ReplacementStrategy::SFIFO: return "SFIFO";
        case lang::ReplacementStrategy::LRU: return "LRU";
        case lang::ReplacementStrategy::HLRU: return "HLRU";
        case lang::ReplacementStrategy::CLIMB: return "CLIMB";
        case lang::ReplacementStrategy::QLRU: return "QLRU";
    }
    throw UnsupportedError("network_writer: unnamed replacement strategy");
}

inline const char* polling_to_json(lang::PollingType p) {
    switch (p) {
        case lang::PollingType::EXHAUSTIVE: return "EXHAUSTIVE";
        case lang::PollingType::GATED: return "GATED";
        case lang::PollingType::KLIMITED: return "KLIMITED";
        case lang::PollingType::DECREMENTING: return "DECREMENTING";
    }
    throw UnsupportedError("network_writer: unnamed polling type");
}

inline const char* impatience_to_json(lang::ImpatienceType i) {
    switch (i) {
        case lang::ImpatienceType::RENEGING: return "RENEGING";
        case lang::ImpatienceType::BALKING: return "BALKING";
        case lang::ImpatienceType::RETRIAL: return "RETRIAL";
        case lang::ImpatienceType::NONE: break;
    }
    throw UnsupportedError("network_writer: a patience block with no impatience type");
}

inline const char* balking_to_json(lang::BalkingStrategy b) {
    switch (b) {
        case lang::BalkingStrategy::QUEUE_LENGTH: return "QUEUE_LENGTH";
        case lang::BalkingStrategy::EXPECTED_WAIT: return "EXPECTED_WAIT";
        case lang::BalkingStrategy::COMBINED: return "COMBINED";
        case lang::BalkingStrategy::NONE: break;
    }
    throw UnsupportedError("network_writer: a balking block with no strategy");
}

/**
 * The `routingStrategies` name of a dispatcher.
 *
 * NOT `lang::routing_to_text`, which is the lower-case name `sn.routing` prints
 * in a struct dump. The wire carries the enum CONSTANT, upper case, matching
 * the JAR's `RoutingStrategy` names one for one -- and the readers compare it
 * exactly, so `rrobin` reloads as an unknown strategy rather than as RROBIN.
 */
inline const char* routing_to_json(lang::RoutingStrategy r) {
    switch (r) {
        case lang::RoutingStrategy::RAND: return "RAND";
        case lang::RoutingStrategy::PROB: return "PROB";
        case lang::RoutingStrategy::RROBIN: return "RROBIN";
        case lang::RoutingStrategy::WRROBIN: return "WRROBIN";
        case lang::RoutingStrategy::JSQ: return "JSQ";
        case lang::RoutingStrategy::FIRING: return "FIRING";
        case lang::RoutingStrategy::SQ: return "SQ";
        case lang::RoutingStrategy::SDR: return "SDR";
        case lang::RoutingStrategy::DISABLED: return "DISABLED";
    }
    throw UnsupportedError("network_writer: unnamed routing strategy");
}

inline const char* hetero_to_json(lang::HeteroSchedPolicy h) {
    switch (h) {
        case lang::HeteroSchedPolicy::ORDER: return "ORDER";
        case lang::HeteroSchedPolicy::ALIS: return "ALIS";
        case lang::HeteroSchedPolicy::ALFS: return "ALFS";
        case lang::HeteroSchedPolicy::FAIRNESS: return "FAIRNESS";
        case lang::HeteroSchedPolicy::FSF: return "FSF";
        case lang::HeteroSchedPolicy::RAIS: return "RAIS";
    }
    throw UnsupportedError("network_writer: unnamed heterogeneous scheduling policy");
}

/**
 * The `scheduling` name of a discipline.
 *
 * UPPER CASE, which is not what `lang::sched_to_text` returns: that is the
 * lower-case spelling `sn.sched` prints in a struct dump, while the wire
 * carries the enum constant (`upper(SchedStrategy.toText(...))` in
 * `linemodel_save.m`, matching `jline.lang.constant.SchedStrategy`). A
 * lower-case `fcfs` reloads as an unknown discipline.
 */
inline std::string sched_to_json(lang::SchedStrategy s) {
    std::string out = lang::sched_to_text(s);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[i])));
    return out;
}

/** A dense matrix as the array of row arrays the readers expect. */
template <class T>
json mat_to_json(const Matrix<T>& M) {
    json rows = json::array();
    for (std::size_t i = 0; i < M.rows(); ++i) {
        json row = json::array();
        for (std::size_t j = 0; j < M.cols(); ++j)
            row.push_back(num_traits<T>::to_double(M(i, j)));
        rows.push_back(row);
    }
    return rows;
}

template <class T>
json vec_to_json(const std::vector<T>& v) {
    json a = json::array();
    for (const T& x : v) a.push_back(num_traits<T>::to_double(x));
    return a;
}

/**
 * Rewrite every non-finite number into the wire's own INFINITY SPELLING.
 *
 * JSON has no infinity literal, so `linemodel_save` -- the reference -- writes
 * an infinite scalar as the STRING "Infinity" / "-Infinity" and a NaN as
 * `null`, a rule it applies to EVERY numeric field rather than a named few;
 * `num_from_json` on the reading side decodes exactly that pair. nlohmann
 * silently dumps an infinity as `null` instead, which the reader then takes for
 * a NaN, so a state row that carries the sentinel -- an open class holds an
 * infinite population, so a Source's `initialState` does the moment the model
 * is initialized -- crossed as a hole rather than as a value. Applied once over
 * the finished document, as the JAR and Python writers apply theirs, so no
 * individual field has to remember the rule.
 */
inline void wire_nonfinite(json& j) {
    if (j.is_object() || j.is_array()) {
        for (json::iterator it = j.begin(); it != j.end(); ++it) wire_nonfinite(*it);
        return;
    }
    if (!j.is_number_float()) return;
    const double v = j.get<double>();
    if (std::isnan(v))
        j = json();
    else if (std::isinf(v))
        j = json(v > 0.0 ? "Infinity" : "-Infinity");
}

/**
 * `x` of a `DiscreteSampler` whose support is the item index, i.e. 1..n.
 *
 * MATLAB `DiscreteSampler(p)` defaults `x = 1:n` and every writer records it;
 * the JAR and Python readers require the key and fail on a file that omits it,
 * so a pmf-only document is unreadable outside this port even though the value
 * is implied.
 */
inline json discrete_support(std::size_t n) {
    json a = json::array();
    for (std::size_t i = 1; i <= n; ++i) a.push_back(static_cast<double>(i));
    return a;
}

/**
 * A distribution as its model.json record, the inverse of `dist_from_json`.
 *
 * The parameter-carrying families are written from `Distrib::params`, which
 * holds the constructor arguments in MATLAB's getParam order -- that is what
 * makes the record reload into the SAME object rather than into a phase-type
 * that shares two moments with it.
 */
template <class T>
json dist_to_json(const lang::Distrib<T>& d) {
    using lang::ProcessType;
    json j;
    j["type"] = lang::process_to_text(d.type);
    const std::vector<T>& p = d.params;
    auto need = [&](std::size_t n, const char* who) {
        if (p.size() < n)
            throw InputError(std::string("network_writer: a ") + who + " carries " +
                             std::to_string(p.size()) + " parameters, and " + std::to_string(n) +
                             " are needed to write it back");
    };
    switch (d.type) {
        case ProcessType::DISABLED:
        case ProcessType::IMMEDIATE:
            return j;  // both are bare tags
        case ProcessType::EXP:
            need(1, "Exp");
            j["params"]["lambda"] = num_traits<T>::to_double(p[0]);
            return j;
        case ProcessType::DET:
            j["params"]["value"] = num_traits<T>::to_double(d.mean);
            return j;
        case ProcessType::ERLANG:
            need(2, "Erlang");
            j["params"]["lambda"] = num_traits<T>::to_double(p[0]);
            j["params"]["k"] = static_cast<long>(std::lround(num_traits<T>::to_double(p[1])));
            return j;
        case ProcessType::HYPEREXP: {
            // The two-branch form stores (p, lambda1, lambda2) and the general
            // one (p_1..p_n, lambda_1..lambda_n); both write the same two
            // vectors, which is the only form on the wire.
            json pv = json::array(), lv = json::array();
            if (p.size() == 3) {
                pv.push_back(num_traits<T>::to_double(p[0]));
                pv.push_back(1.0 - num_traits<T>::to_double(p[0]));
                lv.push_back(num_traits<T>::to_double(p[1]));
                lv.push_back(num_traits<T>::to_double(p[2]));
            } else {
                const std::size_t n = p.size() / 2;
                for (std::size_t i = 0; i < n; ++i) pv.push_back(num_traits<T>::to_double(p[i]));
                for (std::size_t i = 0; i < n; ++i) lv.push_back(num_traits<T>::to_double(p[n + i]));
            }
            j["params"]["p"] = pv;
            j["params"]["lambda"] = lv;
            return j;
        }
        case ProcessType::COXIAN:
        case ProcessType::COX2: {
            const std::size_t n = p.size() / 2;
            json mu = json::array(), phi = json::array();
            for (std::size_t i = 0; i < n; ++i) mu.push_back(num_traits<T>::to_double(p[i]));
            for (std::size_t i = 0; i < n; ++i) phi.push_back(num_traits<T>::to_double(p[n + i]));
            // Coxian is the only name the readers reconstruct from (mu, phi);
            // Cox2 is the two-phase constructor and shares the representation.
            j["type"] = "Coxian";
            j["params"]["mu"] = mu;
            j["params"]["phi"] = phi;
            return j;
        }
        case ProcessType::PH:
        case ProcessType::APH:
        case ProcessType::ME: {
            json alpha = json::array();
            for (const T& v : p) alpha.push_back(num_traits<T>::to_double(v));
            if (d.type == ProcessType::ME) {
                j["params"]["alpha"] = alpha;
                j["params"]["A"] = mat_to_json(d.D0);
            } else {
                j["ph"]["alpha"] = alpha;
                j["ph"]["T"] = mat_to_json(d.D0);
            }
            return j;
        }
        case ProcessType::MAP:
            j["map"]["D0"] = mat_to_json(d.D0);
            j["map"]["D1"] = mat_to_json(d.D1);
            return j;
        case ProcessType::RAP:
            j["params"]["H0"] = mat_to_json(d.D0);
            j["params"]["H1"] = mat_to_json(d.D1);
            return j;
        case ProcessType::DMAP:
            j["params"]["D0"] = mat_to_json(d.D0);
            j["params"]["D1"] = mat_to_json(d.D1);
            return j;
        case ProcessType::MMPP2: {
            // Recovered from the pair rather than from `params`, which
            // `map_dist` does not fill: lambda_i is the arrival rate in phase i
            // and sigma_i the modulating rate out of it.
            if (d.D0.rows() != 2)
                throw InputError("network_writer: an MMPP2 must be of order two");
            j["params"]["lambda0"] = num_traits<T>::to_double(d.D1(0, 0));
            j["params"]["lambda1"] = num_traits<T>::to_double(d.D1(1, 1));
            j["params"]["sigma0"] = num_traits<T>::to_double(d.D0(0, 1));
            j["params"]["sigma1"] = num_traits<T>::to_double(d.D0(1, 0));
            return j;
        }
        case ProcessType::MMAP: {
            j["mmap"]["D0"] = mat_to_json(d.D0);
            json blocks = json::array();
            for (const Matrix<T>& Dk : d.Dmark) blocks.push_back(mat_to_json(Dk));
            j["mmap"]["D1k"] = blocks;
            return j;
        }
        case ProcessType::BMAP: {
            json blocks = json::array();
            blocks.push_back(mat_to_json(d.D0));
            for (const Matrix<T>& Dk : d.Dmark) blocks.push_back(mat_to_json(Dk));
            j["params"]["D"] = blocks;
            return j;
        }
        case ProcessType::UNIFORM:
            need(2, "Uniform");
            j["params"]["a"] = num_traits<T>::to_double(p[0]);
            j["params"]["b"] = num_traits<T>::to_double(p[1]);
            return j;
        case ProcessType::PARETO:
            need(2, "Pareto");
            j["params"]["alpha"] = num_traits<T>::to_double(p[0]);
            j["params"]["scale"] = num_traits<T>::to_double(p[1]);
            return j;
        case ProcessType::GAMMA:
            need(2, "Gamma");
            j["params"]["alpha"] = num_traits<T>::to_double(p[0]);
            j["params"]["beta"] = num_traits<T>::to_double(p[1]);
            return j;
        case ProcessType::WEIBULL:
            // `alpha` is the SCALE here and the SHAPE in Pareto above; that is
            // the reference's own key naming and the one trap in this table.
            need(2, "Weibull");
            j["params"]["alpha"] = num_traits<T>::to_double(p[0]);
            j["params"]["beta"] = num_traits<T>::to_double(p[1]);
            return j;
        case ProcessType::LOGNORMAL:
            need(2, "Lognormal");
            j["params"]["mu"] = num_traits<T>::to_double(p[0]);
            j["params"]["sigma"] = num_traits<T>::to_double(p[1]);
            return j;
        case ProcessType::DUNIFORM:
            need(2, "DiscreteUniform");
            j["params"]["min"] = num_traits<T>::to_double(p[0]);
            j["params"]["max"] = num_traits<T>::to_double(p[1]);
            return j;
        case ProcessType::BERNOULLI:
            need(1, "Bernoulli");
            j["params"]["p"] = num_traits<T>::to_double(p[0]);
            return j;
        case ProcessType::BINOMIAL:
            need(2, "Binomial");
            j["params"]["n"] = static_cast<long>(std::lround(num_traits<T>::to_double(p[0])));
            j["params"]["p"] = num_traits<T>::to_double(p[1]);
            return j;
        case ProcessType::POISSON:
            need(1, "Poisson");
            j["params"]["lambda"] = num_traits<T>::to_double(p[0]);
            return j;
        case ProcessType::GEOMETRIC:
            need(1, "Geometric");
            j["params"]["p"] = num_traits<T>::to_double(p[0]);
            return j;
        case ProcessType::ZIPF:
            need(2, "Zipf");
            j["params"]["s"] = num_traits<T>::to_double(p[0]);
            j["params"]["n"] = static_cast<long>(std::lround(num_traits<T>::to_double(p[1])));
            return j;
        case ProcessType::DISCRETESAMPLER:
            j["params"]["p"] = vec_to_json(p);
            if (!d.trace.empty()) j["params"]["x"] = vec_to_json(d.trace);
            return j;
        case ProcessType::EMPIRICALCDF:
            j["params"]["x"] = vec_to_json(d.trace);
            j["params"]["F"] = vec_to_json(p);
            return j;
        case ProcessType::REPLAYER:
            // The trace itself has no wire form: every writer emits the file
            // PATH, plus the mean beside it as the documented fallback for a
            // reader on another machine where the path does not resolve. A
            // Replayer built from in-memory samples has no path, and then only
            // the moments can be written -- which reloads as a distribution
            // matching them, not as the trace, exactly as the reference warns.
            if (!d.trace_file.empty()) j["params"]["fileName"] = d.trace_file;
            j["params"]["mean"] = num_traits<T>::to_double(d.mean);
            if (d.trace_file.empty())
                j["params"]["scv"] = num_traits<T>::to_double(d.scv);
            return j;
        case ProcessType::NHPP:
        case ProcessType::MAPT:
        case ProcessType::PHT: {
            j["params"]["breakpoints"] = vec_to_json(d.sched_bp);
            j["params"]["cyclic"] = d.sched_cyclic;
            if (d.type == ProcessType::NHPP) {
                // A one-phase schedule: the rate of segment k is D1[k](0,0).
                json rates = json::array();
                for (const Matrix<T>& D1 : d.sched_D1)
                    rates.push_back(num_traits<T>::to_double(D1(0, 0)));
                j["params"]["rates"] = rates;
                return j;
            }
            json A = json::array(), B = json::array();
            for (const Matrix<T>& M : d.sched_D0) A.push_back(mat_to_json(M));
            for (const Matrix<T>& M : d.sched_D1) B.push_back(mat_to_json(M));
            // A PHt is STORED converted to (D0, D1) -- see the Distrib header --
            // so it is written in the MAPt spelling it reloads identically from.
            j["type"] = "MAPt";
            j["params"]["D0"] = A;
            j["params"]["D1"] = B;
            return j;
        }
        case ProcessType::PRIOR:
            throw UnsupportedError(
                "network_writer: a Prior is a set of alternative MODELS rather than one law; save "
                "the design point SolverUQ built, not the design");
        default:
            break;
    }
    throw UnsupportedError(std::string("network_writer: distribution family '") +
                           lang::process_to_text(d.type) + "' has no model.json record");
}

/**
 * The per-class cutoffs a materialized rate table is swept over: a closed
 * class's population, and 10 for an open one.
 *
 * The same rule as `linemodel_save.m:157-163`, and it has to be, or the tables
 * this writer emits would be keyed over a different lattice than the reference's
 * for the same model and the two files would not compare.
 */
template <class T>
std::vector<int> lattice_cutoffs(const qn::NetworkStruct<T>& sn) {
    std::vector<int> cut(sn.classes.size(), 10);
    for (std::size_t r = 0; r < sn.classes.size(); ++r)
        if (sn.classes[r].type == lang::JobClassType::CLOSED &&
            std::isfinite(sn.classes[r].population))
            cut[r] = static_cast<int>(std::lround(sn.classes[r].population));
    return cut;
}

/** Enumerate the box lattice 0 <= n(r) <= cut(r), calling `visit(counts, key)`. */
template <class F>
void for_each_lattice_point(const std::vector<int>& cut, const F& visit) {
    std::size_t total = 1;
    for (int c : cut) total *= static_cast<std::size_t>(c + 1);
    const std::size_t K = cut.size();
    for (std::size_t i = 0; i < total; ++i) {
        std::vector<int> cnt(K, 0);
        std::size_t li = i;
        for (std::size_t d = 0; d < K; ++d) {
            cnt[d] = static_cast<int>(li % static_cast<std::size_t>(cut[d] + 1));
            li /= static_cast<std::size_t>(cut[d] + 1);
        }
        std::string key;
        int tot = 0;
        for (std::size_t d = 0; d < K; ++d) {
            if (d) key += ',';
            key += std::to_string(cnt[d]);
            tot += cnt[d];
        }
        if (tot == 0) continue;  // the empty state is omitted, as the reference omits it
        visit(cnt, key);
    }
}

}  // namespace detail

/**
 * `qn::NetworkStruct` -> the model.json `model` object.
 *
 * Takes the struct rather than the `qn::Network` builder because that is what a
 * solver holds, and because `get_struct()` is what the builder hands out.
 */
template <class T>
detail::json network_to_json(const qn::NetworkStruct<T>& sn) {
    using detail::json;
    typedef qn::NetworkStruct<T> SN;
    const std::size_t K = sn.classes.size();
    const std::vector<int> cut = detail::lattice_cutoffs(sn);

    json model;
    model["type"] = "Network";
    model["name"] = sn.name;
    // The directory every Logger writes into; model-level because that is where
    // the reference keeps it and because a Logger's constructor refuses without it.
    if (!sn.log_path.empty()) model["logPath"] = sn.log_path;

    // -- classes --------------------------------------------------------------
    json classes = json::array();
    for (std::size_t r = 0; r < K; ++r) {
        const qn::JobClass& c = sn.classes[r];
        json cj;
        cj["name"] = c.name;
        // A G-network SIGNAL is written as its own wire type with the kind it
        // removes under: written back as an ordinary class it becomes inert,
        // and the reloaded model is a queueing network with no removals at all.
        const bool sig = r < sn.issignal.size() && sn.issignal[r];
        if (sig) {
            cj["type"] = "Signal";
            cj["openOrClosed"] = c.type == lang::JobClassType::CLOSED ? "Closed" : "Open";
            const lang::SignalType st = sn.signaltype[r];
            cj["signalType"] = st == lang::SignalType::REPLY
                                   ? "reply"
                                   : (st == lang::SignalType::CATASTROPHE ? "catastrophe"
                                                                          : "negative");
            if (sn.signaltarget[r] >= 1 && sn.signaltarget[r] <= K)
                cj["targetClass"] = sn.classes[sn.signaltarget[r] - 1].name;
            const lang::RemovalPolicy rp = sn.signalrempolicy[r];
            if (rp != lang::RemovalPolicy::RANDOM)
                cj["removalPolicy"] = rp == lang::RemovalPolicy::FCFS ? "FCFS" : "LCFS";
            if (!sn.signalremdist[r].empty()) {
                json pv = json::array();
                for (const T& v : sn.signalremdist[r]) pv.push_back(num_traits<T>::to_double(v));
                json xv = json::array();
                // THE SUPPORT STARTS AT ZERO, unlike every other pmf written
                // here: `signalremdist` is indexed BY BATCH SIZE, and a batch of
                // zero is a real outcome. `discrete_support` numbers items from
                // 1, so using it shifted every batch by one job on the way back
                // through any reader that honours `x` -- MATLAB's included.
                for (std::size_t b = 0; b < pv.size(); ++b) xv.push_back(static_cast<double>(b));
                json rd;
                rd["type"] = "DiscreteSampler";
                rd["params"]["p"] = pv;
                rd["params"]["x"] = xv;
                cj["removalDistribution"] = rd;
            }
        } else {
            cj["type"] = c.type == lang::JobClassType::CLOSED
                             ? (c.self_looping ? "SelfLooping" : "Closed")
                             : "Open";
        }
        if (c.type == lang::JobClassType::CLOSED) {
            cj["population"] = c.population;
            if (c.refstat == 0 || c.refstat > sn.station_to_node.size())
                throw InputError("network_writer: class '" + c.name +
                                 "' names no reference station, which a closed class must have");
            cj["refNode"] = sn.nodes[sn.station_to_node[c.refstat - 1] - 1].name;
        }
        if (c.prio != 0) cj["priority"] = c.prio;
        if (c.immfeed) cj["immediateFeedback"] = true;
        if (c.is_ref_class) cj["isReferenceClass"] = true;
        if (std::isfinite(c.deadline)) cj["deadline"] = c.deadline;
        if (c.spawn >= 1 && c.spawn <= K) cj["spawnClass"] = sn.classes[c.spawn - 1].name;
        if (r < sn.syncreply.size() && sn.syncreply[r] >= 1 && sn.syncreply[r] <= K)
            cj["replySignalClass"] = sn.classes[sn.syncreply[r] - 1].name;
        classes.push_back(cj);
    }
    model["classes"] = classes;

    // -- nodes ----------------------------------------------------------------
    json nodes = json::array();
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        const qn::NodeDef& nd = sn.nodes[i];
        const std::size_t ind = i + 1, ist = nd.station;
        json nj;
        nj["name"] = nd.name;
        nj["type"] = detail::node_type_to_json(nd.nodetype);

        if (nd.nodetype == lang::NodeType::Fork && nd.tasks_per_link > 1.0)
            nj["tasksPerLink"] = nd.tasks_per_link;
        // VARIABLE FORKING LEVELS, in the record shape `linemodel_save.m:509-553`
        // writes: one {dest, class, value} per link rather than a matrix, since
        // that is what the reader replays through the setters.
        //
        // Every link a class takes is listed, not only the ones the model
        // overrode. The reader seeds an untouched link from `tasksPerLink`,
        // which is the MEAN once an override exists, so a partial list would
        // reload a fork that emits a count it was never given.
        if (nd.nodetype == lang::NodeType::Fork) {
            const qn::ForkParam<T>* fp = sn.fork_param_of(ind);
            if (fp != 0) {
                json by_dest = json::array(), by_prob = json::array(), by_dist = json::array();
                for (std::size_t k = 0; k < fp->fan_out_link.rows(); ++k)
                    for (std::size_t r = 0; r < fp->fan_out_link.cols() && r < K; ++r) {
                        const double p = num_traits<T>::to_double(fp->fan_out_prob(k, r));
                        if (p == 0.0) continue;  // link this class does not take
                        json rec;
                        rec["dest"] = sn.nodes[k].name;
                        rec["class"] = r + 1;
                        rec["value"] = num_traits<T>::to_double(fp->fan_out_link(k, r));
                        by_dest.push_back(rec);
                        if (p != 1.0) {
                            json pr;
                            pr["dest"] = sn.nodes[k].name;
                            pr["class"] = r + 1;
                            pr["value"] = p;
                            by_prob.push_back(pr);
                        }
                        const lang::Distrib<T>& d = fp->fan_out_dist[k][r];
                        if (!d.disabled) {
                            json dr;
                            dr["dest"] = sn.nodes[k].name;
                            dr["class"] = r + 1;
                            json pv = json::array(), xv = json::array();
                            for (std::size_t e = 0; e < d.params.size(); ++e) {
                                pv.push_back(num_traits<T>::to_double(d.params[e]));
                                xv.push_back(d.trace.empty()
                                                 ? static_cast<double>(e + 1)
                                                 : num_traits<T>::to_double(d.trace[e]));
                            }
                            dr["p"] = pv;
                            dr["x"] = xv;
                            by_dist.push_back(dr);
                        }
                    }
                if (!by_dest.empty()) nj["fanOutByDest"] = by_dest;
                if (!by_prob.empty()) nj["fanOutProb"] = by_prob;
                if (!by_dist.empty()) nj["fanOutDist"] = by_dist;
            }
        }
        // A Logger with no file name exports as a tunnel that logs nowhere, so
        // the name is what makes the node mean anything on reload.
        if (nd.nodetype == lang::NodeType::Logger && !nd.logger.file_name.empty())
            nj["fileName"] = nd.logger.file_name;
        if (nd.nodetype == lang::NodeType::Join) {
            typename std::map<std::size_t, std::pair<std::size_t, std::size_t> >::const_iterator fj;
            for (std::size_t f = 0; f < sn.fj.size(); ++f)
                if (sn.fj[f].second == ind)
                    nj["forkNode"] = sn.nodes[sn.fj[f].first - 1].name;
            typename std::map<std::size_t, typename SN::JoinDecl>::const_iterator jd =
                sn.joindecl.find(ind);
            if (jd != sn.joindecl.end()) {
                if (jd->second.strategy == lang::JoinStrategy::PARTIAL)
                    nj["joinStrategy"] = "PARTIAL";
                if (jd->second.quorum > 0) nj["joinQuorum"] = jd->second.quorum;
            }
        }
        {
            typename std::map<std::size_t, Matrix<T> >::const_iterator cs = sn.csmatrix.find(ind);
            if (cs != sn.csmatrix.end()) {
                json csm;
                for (std::size_t r = 0; r < K; ++r) {
                    json row;
                    for (std::size_t s = 0; s < K; ++s) {
                        const double v = num_traits<T>::to_double(cs->second(r, s));
                        if (v != 0.0) row[sn.classes[s].name] = v;
                    }
                    if (!row.empty()) csm[sn.classes[r].name] = row;
                }
                nj["classSwitchMatrix"] = csm;
            }
        }
        {
            typename std::map<std::size_t, qn::CacheParam<T> >::const_iterator cp =
                sn.nodeparam.find(ind);
            if (cp != sn.nodeparam.end()) {
                const qn::CacheParam<T>& c = cp->second;
                nj["numItems"] = c.nitems;
                nj["itemLevelCap"] = c.itemcap;
                nj["replacementStrategy"] = detail::replacement_to_json(c.replacestrat);
                if (num_traits<T>::to_double(c.qlru) != 1.0)
                    nj["admissionProb"] = num_traits<T>::to_double(c.qlru);
                if (!c.itemsize.empty()) nj["itemSizes"] = c.itemsize;
                if (!c.costcap.empty()) {
                    if (c.costcapglobal) nj["costCaps"] = c.costcap.empty() ? 0 : c.costcap[0];
                    else nj["costCaps"] = c.costcap;
                }
                json pop, hit, miss, itemcls;
                for (std::size_t r = 0; r < K && r < c.pread.size(); ++r) {
                    if (c.pread[r].empty()) continue;
                    json pj;
                    pj["type"] = "DiscreteSampler";
                    pj["params"]["p"] = detail::vec_to_json(c.pread[r]);
                    pj["params"]["x"] = detail::discrete_support(c.pread[r].size());
                    pop[sn.classes[r].name] = pj;
                }
                for (std::size_t r = 0; r < K && r < c.hitclass.size(); ++r)
                    if (c.hitclass[r] != 0)
                        hit[sn.classes[r].name] = sn.classes[c.hitclass[r] - 1].name;
                for (std::size_t r = 0; r < K && r < c.missclass.size(); ++r)
                    if (c.missclass[r] != 0)
                        miss[sn.classes[r].name] = sn.classes[c.missclass[r] - 1].name;
                for (std::size_t r = 0; r < K && r < c.classitem.size(); ++r)
                    if (c.classitem[r] != 0)
                        itemcls[sn.classes[r].name] = static_cast<double>(c.classitem[r]);
                if (!pop.empty()) nj["popularity"] = pop;
                if (!hit.empty()) nj["hitClass"] = hit;
                if (!miss.empty()) nj["missClass"] = miss;
                if (!itemcls.empty()) nj["itemClass"] = itemcls;
                if (!c.accost.empty()) {
                    json ap = json::array();
                    for (const std::vector<Matrix<T> >& per_class : c.accost) {
                        json row = json::array();
                        for (const Matrix<T>& g : per_class)
                            row.push_back(g.rows() == 0 ? json::array() : detail::mat_to_json(g));
                        ap.push_back(row);
                    }
                    nj["accessProb"] = ap;
                }
                if (!c.initstate.empty()) nj["initialState"] = detail::vec_to_json(c.initstate);
                if (c.retrieval_capacity > 0 || !c.retrieval_queues.empty()) {
                    json rs;
                    rs["capacity"] = c.retrieval_capacity;
                    json by;
                    for (const auto& kv : c.retrieval_queues) {
                        json entry;
                        json qs = json::array();
                        for (std::size_t q : kv.second) qs.push_back(sn.nodes[q - 1].name);
                        entry["queues"] = qs;
                        json items;
                        for (std::size_t it = 0; it < c.retrieval_classes.size(); ++it)
                            if (kv.first < c.retrieval_classes[it].size() &&
                                c.retrieval_classes[it][kv.first] != 0)
                                items[std::to_string(it)] =
                                    sn.classes[c.retrieval_classes[it][kv.first] - 1].name;
                        entry["items"] = items;
                        by[sn.classes[kv.first].name] = entry;
                    }
                    rs["byClass"] = by;
                    nj["retrievalSystem"] = rs;
                }
            }
        }
        {
            typename std::map<std::size_t, std::vector<T> >::const_iterator im =
                sn.initmarking.find(ind);
            if (im != sn.initmarking.end()) nj["initialState"] = detail::vec_to_json(im->second);
            typename std::map<std::size_t, std::vector<T> >::const_iterator sp =
                sn.stateprior.find(ind);
            typename std::map<std::size_t, Matrix<T> >::const_iterator ss = sn.statespace.find(ind);
            if (sp != sn.stateprior.end() && ss != sn.statespace.end()) {
                nj["stateSpace"] = detail::mat_to_json(ss->second);
                nj["statePrior"] = detail::vec_to_json(sp->second);
            }
        }

        if (ist != 0) {
            const qn::Station<T>& st = sn.stations[ist - 1];
            if (nd.nodetype == lang::NodeType::Queue || nd.nodetype == lang::NodeType::Delay ||
                nd.nodetype == lang::NodeType::Place)
                nj["scheduling"] = detail::sched_to_json(st.sched);
            if (nd.nodetype == lang::NodeType::Queue && std::isfinite(st.nservers) &&
                st.nservers > 1.0)
                nj["servers"] = st.nservers;
            if (std::isfinite(st.cap) && st.cap > 0) nj["buffer"] = st.cap;

            json svc, cc, dr, sp, imf;
            for (std::size_t r = 0; r < K; ++r) {
                if (ist - 1 < sn.service.size() && r < sn.service[ist - 1].size() &&
                    !sn.service[ist - 1][r].disabled)
                    svc[sn.classes[r].name] = detail::dist_to_json(sn.service[ist - 1][r]);
                if (r < st.classcap.size() && std::isfinite(st.classcap[r]))
                    cc[sn.classes[r].name] = st.classcap[r];
                // 0 is the "the model declared none" sentinel of the station's
                // own vector, not a DropStrategy; `refresh_capacity` reads it
                // as `user_rule` for exactly that reason. Naming it would have
                // to invent a rule, so the entry is simply not written.
                if (r < st.droprule.size() && st.droprule[r] != 0)
                    dr[sn.classes[r].name] =
                        detail::drop_to_json(static_cast<lang::DropStrategy>(st.droprule[r]));
                if (r < st.immfeed.size() && st.immfeed[r]) imf[sn.classes[r].name] = true;
            }
            if (st.sched == lang::SchedStrategy::DPS || st.sched == lang::SchedStrategy::GPS ||
                st.sched == lang::SchedStrategy::DPSPRIO || st.sched == lang::SchedStrategy::GPSPRIO)
                for (std::size_t r = 0; r < K && r < st.schedparam.size(); ++r)
                    sp[sn.classes[r].name] = num_traits<T>::to_double(st.schedparam[r]);
            if (!svc.empty()) nj["service"] = svc;
            if (!cc.empty()) nj["classCap"] = cc;
            if (!dr.empty()) nj["dropRule"] = dr;
            if (!sp.empty()) nj["schedParams"] = sp;
            if (!imf.empty()) nj["immediateFeedback"] = imf;

            if (!st.lldscaling.empty()) {
                json ld;
                ld["type"] = "loadDependent";
                ld["scaling"] = detail::vec_to_json(st.lldscaling);
                nj["loadDependence"] = ld;
            }
            // The two rate-scaling handles, materialized over the box lattice:
            // the only way a function reaches the wire.
            for (int kind = 0; kind < 2; ++kind) {
                const lang::CdScaling<T>& fun = kind == 0 ? st.cdscaling : st.jdscaling;
                if (!static_cast<bool>(fun)) continue;
                json blk, tbl;
                blk["type"] = kind == 0 ? "classDependent" : "jointDependent";
                blk["cutoffs"] = cut;
                detail::for_each_lattice_point(cut, [&](const std::vector<int>& cnt,
                                                        const std::string& key) {
                    std::vector<T> n(K, num_traits<T>::from_int(0));
                    for (std::size_t r = 0; r < K; ++r)
                        n[r] = num_traits<T>::from_int(static_cast<long>(cnt[r]));
                    tbl[key] = detail::vec_to_json(fun(n));
                });
                blk["scaling"] = tbl;
                blk["peak"] = detail::vec_to_json(kind == 0 ? st.cdscalingpeak : st.jdscalingpeak);
                nj[kind == 0 ? "classDependence" : "jointDependence"] = blk;
            }
            {
                typename std::map<std::size_t, typename SN::PasParam>::const_iterator pas =
                    sn.pasparam.find(ist);
                if (pas != sn.pasparam.end() && static_cast<bool>(pas->second.svc_rate_fun)) {
                    json tbl;
                    detail::for_each_lattice_point(
                        cut, [&](const std::vector<int>& cnt, const std::string& key) {
                            std::vector<std::size_t> micro;
                            for (std::size_t r = 0; r < K; ++r)
                                for (int c = 0; c < cnt[r]; ++c) micro.push_back(r + 1);
                            tbl[key] = num_traits<T>::to_double(pas->second.svc_rate_fun(micro));
                        });
                    nj["oiServiceRate"] = tbl;
                    nj["oiCutoffs"] = cut;
                    bool any_swap = false;
                    for (const std::vector<bool>& row : pas->second.swap_graph)
                        for (bool v : row) any_swap = any_swap || v;
                    if (any_swap) {
                        json sg = json::array();
                        for (const std::vector<bool>& row : pas->second.swap_graph) {
                            json rj = json::array();
                            for (bool v : row) rj.push_back(v ? 1 : 0);
                            sg.push_back(rj);
                        }
                        nj["swapGraph"] = sg;
                    }
                }
            }
            if (!st.polling_type.empty()) {
                nj["pollingType"] = detail::polling_to_json(st.polling_type[0]);
                if (st.polling_par != 0) nj["pollingPar"] = st.polling_par;
            }
            {
                json so = json::array();
                for (std::size_t r = 0; r < K && r < st.switchover.size(); ++r) {
                    if (st.switchover[r].disabled) continue;
                    json e;
                    e["from"] = sn.classes[r].name;
                    e["distribution"] = detail::dist_to_json(st.switchover[r]);
                    so.push_back(e);
                }
                if (!so.empty()) nj["switchoverTimes"] = so;
            }
            {
                typename std::map<std::size_t, qn::SetupDelayOffParam<T> >::const_iterator sd =
                    sn.setupparam.find(ist);
                if (sd != sn.setupparam.end()) {
                    json su, doff;
                    for (std::size_t r = 0; r < K && r < sd->second.setup.size(); ++r) {
                        if (sd->second.setup[r].disabled) continue;
                        su[sn.classes[r].name] = detail::dist_to_json(sd->second.setup[r]);
                        doff[sn.classes[r].name] = detail::dist_to_json(sd->second.delayoff[r]);
                    }
                    if (!su.empty()) {
                        nj["setupTime"] = su;
                        nj["delayOffTime"] = doff;
                    }
                }
            }
            {
                // Server breakdown, in `linemodel_save`'s shape: the two clocks
                // always, and `downService` keyed by class name only where a
                // degraded rate was declared. A zero rate is an ABSENT entry and
                // not a zero-mean distribution, which has no reading.
                typename std::map<std::size_t, qn::BreakdownParam<T> >::const_iterator bd =
                    sn.breakdownparam.find(ist);
                if (bd != sn.breakdownparam.end()) {
                    json bj;
                    bj["failure"] = detail::dist_to_json(bd->second.failure);
                    bj["repair"] = detail::dist_to_json(bd->second.repair);
                    json ds;
                    for (std::size_t r = 0; r < K && r < bd->second.down_service_rates.size(); ++r) {
                        const double rate =
                            num_traits<T>::to_double(bd->second.down_service_rates[r]);
                        if (!(rate > 0.0)) continue;
                        ds[sn.classes[r].name] = detail::dist_to_json(
                            lang::Distrib<T>::exp_rate(rate));
                    }
                    if (!ds.empty()) bj["downService"] = ds;
                    nj["breakdown"] = bj;
                }
            }
            {
                typename std::map<std::size_t, qn::RetrialParam<T> >::const_iterator rt =
                    sn.retrialparam.find(ist);
                if (rt != sn.retrialparam.end()) {
                    json rj;
                    for (std::size_t r = 0; r < K && r < rt->second.retrial_proc.size(); ++r) {
                        if (rt->second.retrial_proc[r].disabled) continue;
                        json e;
                        e["delay"] = detail::dist_to_json(rt->second.retrial_proc[r]);
                        e["maxAttempts"] = rt->second.max_attempts[r];
                        rj[sn.classes[r].name] = e;
                    }
                    if (!rj.empty()) nj["retrial"] = rj;
                }
            }
            {
                json pat, orb, brp, blk;
                for (std::size_t r = 0; r < K; ++r) {
                    if (r < st.patience.size() && !st.patience[r].disabled) {
                        json e;
                        e["distribution"] = detail::dist_to_json(st.patience[r]);
                        if (r < st.impatience.size() &&
                            st.impatience[r] != lang::ImpatienceType::NONE)
                            e["impatienceType"] = detail::impatience_to_json(st.impatience[r]);
                        pat[sn.classes[r].name] = e;
                    }
                    if (r < st.orbit_impatience.size() && !st.orbit_impatience[r].disabled)
                        orb[sn.classes[r].name] = detail::dist_to_json(st.orbit_impatience[r]);
                    if (r < st.batch_reject.size() &&
                        num_traits<T>::to_double(st.batch_reject[r]) > 0)
                        brp[sn.classes[r].name] = num_traits<T>::to_double(st.batch_reject[r]);
                    if (r < st.balking.size() &&
                        st.balking[r].strategy != lang::BalkingStrategy::NONE) {
                        json e;
                        e["strategy"] = detail::balking_to_json(st.balking[r].strategy);
                        json ths = json::array();
                        for (const typename qn::Station<T>::BalkingThreshold& th :
                             st.balking[r].thresholds) {
                            json tj;
                            tj["minJobs"] = th.min_jobs;
                            tj["maxJobs"] = th.max_jobs;
                            tj["probability"] = num_traits<T>::to_double(th.probability);
                            ths.push_back(tj);
                        }
                        e["thresholds"] = ths;
                        blk[sn.classes[r].name] = e;
                    }
                }
                if (!pat.empty()) nj["patience"] = pat;
                if (!orb.empty()) nj["orbitImpatience"] = orb;
                if (!brp.empty()) nj["batchRejectProb"] = brp;
                if (!blk.empty()) nj["balking"] = blk;
            }
            if (!st.server_types.empty()) {
                json sts = json::array();
                for (const typename qn::Station<T>::ServerType& t : st.server_types) {
                    json tj;
                    tj["name"] = t.name;
                    tj["count"] = t.count;
                    if (!t.compatible.empty()) {
                        json cn = json::array();
                        for (std::size_t r = 0; r < K && r < t.compatible.size(); ++r)
                            if (t.compatible[r]) cn.push_back(sn.classes[r].name);
                        tj["compatibleClasses"] = cn;
                    }
                    json sv;
                    for (std::size_t r = 0; r < K && r < t.service.size(); ++r)
                        if (!t.service[r].disabled)
                            sv[sn.classes[r].name] = detail::dist_to_json(t.service[r]);
                    if (!sv.empty()) tj["service"] = sv;
                    sts.push_back(tj);
                }
                nj["serverTypes"] = sts;
                if (st.hetero_policy != lang::HeteroSchedPolicy::ORDER)
                    nj["heteroSchedPolicy"] = detail::hetero_to_json(st.hetero_policy);
            }
            {
                json par;
                for (std::size_t r = 0; r < K && r < st.server_parallelism.size(); ++r)
                    if (st.server_parallelism[r] > 1)
                        par[sn.classes[r].name] = st.server_parallelism[r];
                if (!par.empty()) nj["serverParallelism"] = par;
            }
            {
                json ab;
                for (std::size_t r = 0; r < K && r < st.arrival_batch.size(); ++r)
                    if (!st.arrival_batch[r].disabled)
                        ab[sn.classes[r].name] = detail::dist_to_json(st.arrival_batch[r]);
                if (!ab.empty()) nj["arrivalBatch"] = ab;
                if (!st.marked_classes.empty()) {
                    json mc = json::array();
                    for (std::size_t c : st.marked_classes) mc.push_back(sn.classes[c - 1].name);
                    nj["markedClasses"] = mc;
                }
                json dd;
                for (std::size_t r = 0; r < K && r < st.departure_discipline.size(); ++r)
                    if (st.departure_discipline[r] != lang::DepartureDiscipline::NORMAL)
                        dd[sn.classes[r].name] = "FIFO";
                if (!dd.empty()) nj["departureDiscipline"] = dd;
            }
        }

        {
            typename std::map<std::size_t, qn::TransitionParam<T> >::const_iterator tp =
                sn.transparam.find(ind);
            if (tp != sn.transparam.end()) {
                json modes = json::array();
                for (std::size_t m = 0; m < tp->second.nmodes; ++m) {
                    json mj;
                    mj["name"] = tp->second.modenames[m];
                    mj["timingStrategy"] =
                        tp->second.timing[m] == lang::TimingStrategy::IMMEDIATE ? "IMMEDIATE"
                                                                                : "TIMED";
                    if (!tp->second.firingproc[m].disabled)
                        mj["distribution"] = detail::dist_to_json(tp->second.firingproc[m]);
                    mj["numServers"] = tp->second.nmodeservers[m];
                    mj["firingPriority"] = tp->second.firingprio[m];
                    mj["firingWeight"] = num_traits<T>::to_double(tp->second.fireweight[m]);
                    const char* kArcKey[3] = {"enablingConditions", "inhibitingConditions",
                                             "firingOutcomes"};
                    for (int which = 0; which < 3; ++which) {
                        const Matrix<T>& row = which == 0   ? tp->second.enabling[m]
                                               : which == 1 ? tp->second.inhibiting[m]
                                                            : tp->second.firing[m];
                        json arcs = json::array();
                        for (std::size_t q = 0; q < row.rows(); ++q)
                            for (std::size_t r = 0; r < row.cols(); ++r) {
                                const double v = num_traits<T>::to_double(row(q, r));
                                // An inhibiting threshold of Inf is "never
                                // blocks" and an enabling or firing count of 0
                                // is no arc.
                                if (which == 1 ? !std::isfinite(v) : v == 0.0) continue;
                                json a;
                                a["node"] = sn.nodes[q].name;
                                // REQUIRED by the reference readers: both index
                                // the arc by (node, CLASS) -- `linemodel_load.m:
                                // 852` reads `ec.class` unguarded and
                                // `linemodel_io.py` does `ic["class"]` -- so an
                                // arc written without it makes the whole
                                // document unloadable there.
                                a["class"] = sn.classes[r].name;
                                a["count"] = v;
                                arcs.push_back(a);
                            }
                        if (!arcs.empty()) mj[kArcKey[which]] = arcs;
                    }
                    // A marking-dependent firing rate is a CLOSURE, so it is
                    // written the only way it can cross JSON: materialized over
                    // the box lattice of the enabling slots, exactly as
                    // `firingdep_scaling_table` in linemodel_save.m does. Dropped
                    // instead, the net reloads at its nominal rate -- a different
                    // marking process, reported without a diagnostic.
                    if (m < tp->second.firingdep.size() && tp->second.firingdep[m]) {
                        const Matrix<T>& enab = tp->second.enabling[m];
                        json slots = json::array();
                        std::vector<std::size_t> slot_node;
                        std::vector<long> caps;
                        // ONE SLOT PER PLACE, not per (place, class): the
                        // closure's own domain is the node-indexed marking that
                        // `state_events.h` hands it, and the reader resolves a
                        // slot by node. The class is the first the mode reads
                        // there, which is the only one a single-class arc set
                        // can name and is what labels the slot for the
                        // reference readers.
                        for (std::size_t q = 0; q < enab.rows(); ++q) {
                            std::size_t rq = enab.cols();
                            for (std::size_t r = 0; r < enab.cols(); ++r)
                                if (num_traits<T>::to_double(enab(q, r)) > 0.0) { rq = r; break; }
                            if (rq == enab.cols()) continue;
                            json sm;
                            sm["node"] = sn.nodes[q].name;
                            sm["class"] = sn.classes[rq].name;
                            slots.push_back(sm);
                            slot_node.push_back(q);
                            // The open-place saturation cutoff of the three
                            // reference writers: an unbounded place would make
                            // the lattice infinite, so it is tabulated to 10.
                            const std::size_t sq = sn.nodes[q].station;
                            const double pc = sq == 0 ? std::numeric_limits<double>::infinity()
                                                      : sn.stations[sq - 1].cap;
                            caps.push_back(std::isfinite(pc) ? std::lround(pc) : 10L);
                        }
                        if (!slots.empty()) {
                            json frm;
                            frm["slots"] = slots;
                            frm["cutoffs"] = caps;
                            json scaling;
                            std::size_t total = 1;
                            for (std::size_t s = 0; s < caps.size(); ++s)
                                total *= static_cast<std::size_t>(caps[s]) + 1;
                            std::vector<T> mk(sn.nodes.size(), num_traits<T>::from_int(0));
                            for (std::size_t li = 0; li < total; ++li) {
                                std::size_t rem = li;
                                std::string key;
                                for (std::size_t s = 0; s < caps.size(); ++s) {
                                    const std::size_t shp = static_cast<std::size_t>(caps[s]) + 1;
                                    const std::size_t c = rem % shp;
                                    rem /= shp;
                                    mk[slot_node[s]] = num_traits<T>::from_int(
                                        static_cast<long>(c));
                                    if (s) key += ',';
                                    key += std::to_string(c);
                                }
                                const double v =
                                    num_traits<T>::to_double(tp->second.firingdep[m](mk));
                                scaling[key] = std::isfinite(v) ? v : 0.0;
                            }
                            for (std::size_t s = 0; s < slot_node.size(); ++s)
                                mk[slot_node[s]] = num_traits<T>::from_int(0);
                            frm["scaling"] = scaling;
                            mj["firingRateDependence"] = frm;
                        }
                    }
                    modes.push_back(mj);
                }
                nj["modes"] = modes;
            }
        }
        nodes.push_back(nj);
    }
    model["nodes"] = nodes;

    // -- routing --------------------------------------------------------------
    json matrix;
    for (const auto& kv : sn.P) {
        const std::size_t r = kv.first.first, s = kv.first.second;
        json from_to;
        for (std::size_t a = 0; a < kv.second.rows(); ++a) {
            json row;
            for (std::size_t b = 0; b < kv.second.cols(); ++b) {
                const double p = num_traits<T>::to_double(kv.second(a, b));
                if (p != 0.0) row[sn.nodes[b].name] = p;
            }
            if (!row.empty()) from_to[sn.nodes[a].name] = row;
        }
        if (!from_to.empty())
            matrix[sn.classes[r - 1].name + "," + sn.classes[s - 1].name] = from_to;
    }
    json routing;
    routing["type"] = "matrix";
    routing["matrix"] = matrix;
    model["routing"] = routing;

    json strategies, weights, params;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        const qn::NodeDef& nd = sn.nodes[i];
        json per_class, per_class_w, per_class_p;
        for (std::size_t r = 0; r < K && r < nd.routing.size(); ++r) {
            // PROB is the matrix itself; RAND is DECLARED not derived, and dropping it wrote JMT's Empirical strategy where the model asks for Random.
            if (nd.routing[r] != lang::RoutingStrategy::PROB)
                per_class[sn.classes[r].name] = detail::routing_to_json(nd.routing[r]);
            if (r < nd.routing_weights.size() && !nd.routing_weights[r].empty()) {
                json dw;
                for (const auto& kv : nd.routing_weights[r]) dw[sn.nodes[kv.first - 1].name] = kv.second;
                per_class_w[sn.classes[r].name] = dw;
            }
            if (r < nd.routing_param.size() && nd.routing_param[r] != 0) {
                json pj;
                pj["d"] = nd.routing_param[r];
                per_class_p[sn.classes[r].name] = pj;
            }
        }
        if (!per_class.empty()) strategies[nd.name] = per_class;
        if (!per_class_w.empty()) weights[nd.name] = per_class_w;
        if (!per_class_p.empty()) params[nd.name] = per_class_p;
    }
    if (!strategies.empty()) model["routingStrategies"] = strategies;
    if (!weights.empty()) model["routingWeights"] = weights;
    if (!params.empty()) model["routingParams"] = params;

    // Krzesinski state-dependent routing. Every centre travels by NODE NAME, so
    // the block is language independent and a node reordering on either side
    // cannot shift a centre. Branch index 1 denotes the complement M-V and is
    // written as an empty list, keeping the paper's own numbering.
    // See _kb/16-state-dependent-routing.md
    if (!sn.sdr_nodes.empty()) {
        const pfqn::SdrStruct& sd = sn.sdr_nodes;  // 0-based NODE indices
        std::size_t cls = 0;
        for (std::size_t r = 0; r < K && cls == 0; ++r)
            if (sn.nodes[sd.entry].routing.size() > r &&
                sn.nodes[sd.entry].routing[r] == lang::RoutingStrategy::SDR)
                cls = r + 1;
        if (cls == 0)
            throw InputError("network_writer: the model declares state-dependent routing at node '" +
                             sn.nodes[sd.entry].name + "' but no class routes SDR there");
        json sdr;
        sdr["entry"] = sn.nodes[sd.entry].name;
        sdr["departure"] = sn.nodes[sd.departure].name;
        sdr["class"] = sn.classes[cls - 1].name;
        json branches = json::array();
        for (std::size_t b = 0; b < sd.branch.size(); ++b) {
            json bn = json::array();
            for (std::size_t q = 0; q < sd.branch[b].size(); ++q)
                bn.push_back(sn.nodes[sd.branch[b][q]].name);
            branches.push_back(bn);
        }
        sdr["branches"] = branches;
        json level = json::array(), Cs = json::array(), dm = json::array();
        for (std::size_t b = 0; b < sd.level.size(); ++b)
            level.push_back(static_cast<double>(sd.level[b]));
        for (std::size_t t = 0; t < sd.C.size(); ++t) Cs.push_back(sd.C[t]);
        for (std::size_t t = 0; t < sd.d.rows(); ++t) {
            json row = json::array();
            for (std::size_t b = 0; b < sd.d.cols(); ++b) row.push_back(sd.d(t, b));
            dm.push_back(row);
        }
        sdr["level"] = level;
        sdr["C"] = Cs;
        sdr["d"] = dm;
        model["stateDepRouting"] = sdr;
    }

    // -- finite capacity regions ---------------------------------------------
    if (!sn.regions.empty()) {
        json fcr = json::array();
        for (std::size_t g = 0; g < sn.regions.size(); ++g) {
            const typename SN::Region& rg = sn.regions[g];
            json rj;
            rj["name"] = rg.name.empty() ? "Region" + std::to_string(g + 1) : rg.name;
            json stations = json::array();
            double global = -1.0, globalmem = -1.0;
            json class_cap, class_size, class_weight;
            for (std::size_t m = 0; m < rg.members.size(); ++m) {
                if (!rg.members[m]) continue;
                json sj;
                sj["node"] = sn.nodes[sn.station_to_node[m] - 1].name;
                stations.push_back(sj);
                global = rg.cap[m][K];
                // The region's global MEMORY budget, which is enforced beside the
                // job cap by the CTMC, the loss-network NC arm and the JMT export
                // rather than folded into it. Omitted, the region reloads with
                // unbounded memory and every one of those three drops a
                // constraint the model declared.
                if (m < rg.maxmem.size()) globalmem = rg.maxmem[m];
                for (std::size_t r = 0; r < K; ++r)
                    if (rg.cap[m][r] != -1.0) class_cap[sn.classes[r].name] = rg.cap[m][r];
            }
            rj["stations"] = stations;
            if (global != -1.0) rj["globalMaxJobs"] = global;
            if (globalmem != -1.0) rj["globalMaxMemory"] = globalmem;
            if (!class_cap.empty()) rj["classMaxJobs"] = class_cap;
            json rule;
            for (std::size_t r = 0; r < K && r < rg.rule.size(); ++r)
                rule[sn.classes[r].name] = detail::drop_to_json(rg.rule[r]);
            rj["dropRule"] = rule;
            for (std::size_t r = 0; r < K; ++r) {
                if (r < rg.size.size() && num_traits<T>::to_double(rg.size[r]) != 1.0)
                    class_size[sn.classes[r].name] = num_traits<T>::to_double(rg.size[r]);
                if (r < rg.weight.size() && num_traits<T>::to_double(rg.weight[r]) != 1.0)
                    class_weight[sn.classes[r].name] = num_traits<T>::to_double(rg.weight[r]);
            }
            // classSize and classWeight are per-STATION on the wire, and the
            // struct holds one vector per region, so they are written onto every
            // member station -- which is the same model the reader rebuilds.
            if (!class_size.empty() || !class_weight.empty())
                for (json& sj : rj["stations"]) {
                    if (!class_size.empty()) sj["classSize"] = class_size;
                    if (!class_weight.empty()) sj["classWeight"] = class_weight;
                }
            if (rg.lincon_A.rows() > 0) {
                rj["constraintA"] = detail::mat_to_json(rg.lincon_A);
                rj["constraintB"] = detail::vec_to_json(rg.lincon_b);
            }
            fcr.push_back(rj);
        }
        model["finiteCapacityRegions"] = fcr;
    }

    // -- global (Whittle) dependence -----------------------------------------
    //
    // phi(n) reads the FULL population matrix, so unlike the per-station
    // classDependence/jointDependence blocks it is materialized over the lattice
    // of the WHOLE network state. Only the (station,class) SLOTS a class can
    // actually occupy carry a coordinate: a Source holds no jobs and a class with
    // zero per-class capacity at a station never appears there. That restriction
    // is lossless, since no DEP or PHASE event ever fires at such a slot.
    if (static_cast<bool>(sn.gdscaling)) {
        const std::size_t M = sn.nstations, K = sn.nclasses;
        const std::vector<double> njobs = sn.njobs();
        const int wcut = sn.gdscalingcutoff;
        std::vector<std::size_t> slot_st, slot_cl;
        std::vector<int> cuts;
        for (std::size_t i = 0; i < M; ++i) {
            if (sn.nodes[sn.station_to_node[i] - 1].nodetype == lang::NodeType::Source) continue;
            for (std::size_t r = 0; r < K; ++r) {
                const double cap = sn.classcap[i][r];
                if (!(cap > 0)) continue;
                int c = std::isfinite(njobs[r]) ? static_cast<int>(std::lround(njobs[r])) : wcut;
                if (std::isfinite(cap)) c = std::min(c, static_cast<int>(std::lround(cap)));
                slot_st.push_back(i);
                slot_cl.push_back(r);
                cuts.push_back(c > 0 ? c : 0);
            }
        }
        const std::size_t P = cuts.size();
        std::size_t total = 1;
        for (std::size_t d = 0; d < P; ++d) {
            total *= static_cast<std::size_t>(cuts[d] + 1);
            if (total > 200000u)
                throw InputError(
                    "the global dependence lattice exceeds the wire limit of 200000 points; lower "
                    "the wireCutoff argument of set_global_dependence, or solve the model "
                    "natively");
        }

        json blk, slots, tbl;
        blk["type"] = "globalDependent";
        std::vector<std::string> station_names(M), class_names(K);
        for (std::size_t i = 0; i < M; ++i)
            station_names[i] = sn.nodes[sn.station_to_node[i] - 1].name;
        for (std::size_t r = 0; r < K; ++r) class_names[r] = sn.classes[r].name;
        blk["stations"] = station_names;
        blk["classes"] = class_names;
        slots = json::array();
        for (std::size_t d = 0; d < P; ++d) {
            json sm;
            sm["station"] = station_names[slot_st[d]];
            sm["class"] = class_names[slot_cl[d]];
            slots.push_back(sm);
        }
        blk["slots"] = slots;
        blk["cutoffs"] = cuts;
        blk["cutoff"] = wcut;

        for (std::size_t li = 0; li < total; ++li) {
            std::size_t rem = li;
            std::vector<int> cnt(P, 0);
            for (std::size_t d = 0; d < P; ++d) {
                cnt[d] = static_cast<int>(rem % static_cast<std::size_t>(cuts[d] + 1));
                rem /= static_cast<std::size_t>(cuts[d] + 1);
            }
            std::vector<T> n(M * K, num_traits<T>::from_int(0));
            for (std::size_t d = 0; d < P; ++d)
                n[slot_st[d] * K + slot_cl[d]] = num_traits<T>::from_int(cnt[d]);
            const std::vector<T> v = sn.gdscaling(n);
            std::vector<double> flat(M * K, 1.0);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    double x;
                    if (v.size() == 1) x = num_traits<T>::to_double(v[0]);
                    else if (v.size() == M) x = num_traits<T>::to_double(v[i]);
                    else x = num_traits<T>::to_double(v[i * K + r]);
                    flat[i * K + r] = std::isfinite(x) ? x : 0.0;
                }
            std::string key;
            if (P == 0) key = "0";
            else
                for (std::size_t d = 0; d < P; ++d) {
                    if (d) key += ',';
                    key += std::to_string(cnt[d]);
                }
            tbl[key] = flat;
        }
        blk["scaling"] = tbl;
        std::vector<double> pk(M * K, 1.0);
        for (std::size_t j = 0; j < pk.size() && j < sn.gdscalingpeak.size(); ++j)
            pk[j] = num_traits<T>::to_double(sn.gdscalingpeak[j]);
        blk["peak"] = pk;
        model["globalDependence"] = blk;
    }

    // -- rewards --------------------------------------------------------------
    //
    // The DECLARATIVE ones only. A reward built from a lambda carries no `kind`
    // and is dropped here, exactly as `linemodel_save` drops it: any record
    // written for it would reload as a different function.
    {
        // IN NAME ORDER, which is the cross-codebase convention and not a
        // cosmetic choice: `rewards2json` sorts for it explicitly, because the
        // JAR holds the rewards in a HashMap and has no insertion order to
        // preserve. Emitting them in struct order made the same model produce a
        // different document here than in the other three, which the JSON arm of
        // the parity harness reports as a difference in the model.
        std::map<std::string, json> by_name;
        for (const typename SN::Reward& rw : sn.reward) {
            if (rw.kind.empty() || rw.node == 0) continue;
            json rj;
            rj["name"] = rw.name;
            rj["type"] = rw.kind;
            rj["node"] = sn.nodes[rw.node - 1].name;
            if (rw.cls != 0) rj["class"] = sn.classes[rw.cls - 1].name;
            by_name[rw.name] = rj;
        }
        json rewards = json::array();
        for (const std::pair<const std::string, json>& kv : by_name) rewards.push_back(kv.second);
        if (!rewards.empty()) model["rewards"] = rewards;
    }
    return model;
}

/** The complete model.json envelope: `{format, version, model}`. */
template <class T>
detail::json network_json_envelope(const qn::NetworkStruct<T>& sn) {
    detail::json root;
    root["format"] = "line-model";
    root["version"] = "1.0";
    root["model"] = network_to_json(sn);
    detail::wire_nonfinite(root);
    return root;
}

/** Write a model.json file, indented as the reference writers indent it. */
template <class T>
void write_network_json(const qn::NetworkStruct<T>& sn, const std::string& path) {
    std::ofstream out(path.c_str());
    if (!out) throw InputError("network_writer: cannot open " + path + " for writing");
    out << network_json_envelope(sn).dump(2) << "\n";
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_NETWORK_WRITER_H
