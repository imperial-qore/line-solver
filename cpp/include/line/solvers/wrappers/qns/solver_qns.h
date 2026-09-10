/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_WRAPPERS_QNS_SOLVER_QNS_H
#define LINE_SOLVERS_WRAPPERS_QNS_SOLVER_QNS_H

/**
 * Port of `@@SolverQNS`, the wrapper around `qnsolver` of the RADS/LQNS
 * distribution.
 *
 * WHY A WRAPPER IS PORTED AT ALL. Every other solver in `cpp/` computes its own
 * answer; this one marshals the model to an external binary and reads the
 * numbers back. It earns its place for the same reason it does in the other
 * three codebases: `qnsolver` is an INDEPENDENT implementation of the multiserver
 * AMVA lineage that `solver_mva` also implements, so it is the cross-check that
 * catches an error common to the port and its reference. It is the only external
 * tool in `cpp/`: LINE ships no copy of it, and every path here refuses by name
 * when the binary is absent rather than answering natively under the QNS label.
 *
 * WHAT THE PORT COVERS. The reference dispatches two ways (`runAnalyzer.m`):
 *
 *   - product-form, or any model with open classes -> marshal to JMVA, run
 *     `qnsolver`, parse and de-aggregate. THIS IS PORTED, in full.
 *   - non-product-form closed -> convert with `QN2LQN` and delegate to
 *     `SolverLQNS`. This is the same layered path as MATLAB and the JAR; LINE
 *     still ships no LQNS binary, so its ordinary availability diagnostic is
 *     preserved when the external tool is absent.
 *
 * METHODS. `conway`, `reiser`, `rolia` and `zhou` are what `qnsolver -m` accepts,
 * and they only take effect on a model that HAS a multiserver station -- the
 * reference passes no `-m` otherwise, and so does this. `suri` and `schmidt` are
 * listed by the reference's `listValidMethods` and reach the tool through the
 * LQNS branch on a non-product-form closed model.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "line/io/jmva_writer.h"
#include "line/io/qn2lqn.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/wrappers/lqns/solver_lqns.h"
#include "line/util/error.h"
#include "line/util/subprocess.h"
#include "line/util/tempdir.h"

namespace line {
namespace qns {

/** `SolverQNS.defaultOptions` plus the two knobs the JMVA document carries. */
struct QnsOptions {
    std::string method = "default";
    /**
     * `options.config.multiserver`. Kept separate from `method` because the
     * reference lets either name the approximation: `method` sets it, and a
     * caller may set it directly while leaving the method at 'default'.
     */
    std::string multiserver;
    std::size_t samples = 10000;  ///< `options.samples`, the JMVA maxSamples
    /** Seconds before a hung `qnsolver` is killed; not positive waits forever. */
    int timeout = 0;
    /** `options.keep`: leave the scratch directory behind, to inspect what was sent. */
    bool keep = false;
};

/** Port of `SolverQNS.listValidMethods`. */
inline std::vector<std::string> list_valid_methods() {
    return {"default", "conway", "rolia", "zhou", "suri", "reiser", "schmidt"};
}

/** The approximations `qnsolver -m` accepts; the rest reach it only via LQNS. */
inline bool is_qnsolver_multiserver(const std::string& m) {
    return m == "conway" || m == "reiser" || m == "rolia" || m == "zhou";
}

/**
 * `SolverQNS.supportsModelMethod`'s structural rules, as the REASON they refuse,
 * empty when the pair is served.
 *
 * A STRING RATHER THAN A THROW, the shape `ba::method_refusal` and
 * `ag::runner_detail::method_refusal` already carry: the AUTO report has to ASK
 * the question without raising, so that a pair it offers is a pair the run
 * accepts. `solver_qns_run` keeps the throws, which are the same two rules
 * stated where the run reaches them.
 *
 * Mirrors matlab/src/solvers/wrappers/QNS/qns_immfeed_refusal.m and
 * qns_multiserver_refusal.m, and the JAR's qnsImmfeedRefusal /
 * qnsMultiserverRefusal.
 */
template <class T>
std::string method_refusal(const qn::NetworkStruct<T>& L, const std::string& method) {
    if (L.has_immediate_feedback())
        return "SolverQNS does not support immediate feedback (sn.immfeed): neither the JMVA "
               "document qnsolver reads nor the LQN QN2LQN writes can keep a self-looping job on "
               "its server. Use SolverCTMC or SolverSSA, whose state space carries the self-loop.";

    // THE RULE IS INSIDE THE MULTISERVER BRANCH: without a multiserver station
    // no -m is emitted and every method name is served by the plain invocation.
    bool has_multiserver = false;
    for (std::size_t i = 0; i < L.nstations && !has_multiserver; ++i) {
        const double c = L.stations[i].nservers;
        if (std::isfinite(c) && c > 1.0) has_multiserver = true;
    }
    if (!has_multiserver) return std::string();
    std::string ms = method.empty() ? std::string("default") : method;
    std::transform(ms.begin(), ms.end(), ms.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (ms == "default" || is_qnsolver_multiserver(ms)) return std::string();
    return "SolverQNS: the multiserver approximation '" + ms +
           "' is one LQNS offers and qnsolver does not: 'qnsolver -m' accepts conway, reiser, "
           "rolia and zhou only; suri and schmidt are available only on the non-product-form "
           "closed SolverLQNS branch.";
}

namespace detail {

/**
 * The argv prefix every invocation uses.
 *
 * LD_LIBRARY_PATH IS STRIPPED, as `solver_qns.m` does on unix: the binaries of
 * the LQNS distribution are linked against the system libstdc++, and a
 * LD_LIBRARY_PATH inherited from a host that ships its own (MATLAB does) makes
 * them fail to load on a GLIBCXX version symbol. Doing it through `env -u`
 * rather than by editing this process's environment keeps the change confined to
 * the child.
 */
inline std::vector<std::string> qnsolver_argv() {
    return {"env", "-u", "LD_LIBRARY_PATH", "qnsolver"};
}

}  // namespace detail

/**
 * Port of `SolverQNS.isAvailable`: a native `qnsolver` binary on the PATH.
 *
 * There is no container fallback, here as in the other three codebases: qnsolver
 * ships with LQNS, whose licence forbids providing the software to another
 * party, so no LINE code path resolves or runs an image carrying it.
 */
inline bool is_available() {
    std::vector<std::string> argv = detail::qnsolver_argv();
    argv.push_back("--help");
    const util::ProcResult r = util::capture(argv, 30);
    return r.exitCode == 0;
}

/** Port of `runAnalyzer`'s method gate. */
inline void check_method(const std::string& method) {
    const std::vector<std::string> valid = list_valid_methods();
    if (std::find(valid.begin(), valid.end(), method) == valid.end())
        throw InputError("SolverQNS: unknown method '" + method +
                         "'; valid methods are default, conway, rolia, zhou, suri, reiser, "
                         "schmidt");
}

/**
 * Port of the method -> `options.config.multiserver` map of `runAnalyzer.m`.
 *
 * `default` resolves to `rolia`, as `runAnalyzer.m` does. THE REFERENCE HAS TWO
 * ENTRY POINTS AND THEY DISAGREE on an explicit config: `runAnalyzer.m`
 * overwrites `options.config.multiserver` unconditionally at method `default`,
 * so a MATLAB caller who sets it through the solver constructor is silently
 * ignored, while `solver_qns.m` called directly honours it and maps its own
 * `default` to CONWAY. This resolves the pair the only way that is a superset of
 * neither: the method wins when named, an explicit config is honoured (the
 * low-level reference's behaviour, and the only way a caller can express the
 * choice at all), and an unset one gives rolia (the high-level reference's).
 * Native Python currently leaves `default` at CONWAY here, which is a genuine
 * divergence from MATLAB and the JAR; see the report accompanying this audit.
 */
inline std::string resolve_multiserver(const QnsOptions& opt) {
    if (opt.method != "default") return opt.method;
    if (!opt.multiserver.empty()) return opt.multiserver;
    return "rolia";
}

namespace detail {

/** One parsed data row of the `qnsolver` CSV output. */
struct ParsedRow {
    std::string station;
    std::vector<double> Q, W, U, Tp;
};

inline double parse_double_or_zero(const std::string& s) {
    if (s.empty()) return 0.0;
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str()) return 0.0;
    return v;
}

/**
 * Split a `qnsolver` result row into its four metric blocks.
 *
 * THE STRIDE IS DECIDED BY THE CHAIN COUNT, not the class count. `qnsolver`
 * writes a per-chain column followed by an aggregate one for each of Q, R, U and
 * X when there is more than one chain, and drops the aggregate entirely at one
 * chain -- confirmed against the binary: one chain gives
 * `Station, $Q, $R, $U, $X` and two give `$Q(Chain01), $Q(Chain02), $Q, ...`.
 * MATLAB and the JAR both switch on `sn.nclasses == 1`, which agrees only while
 * classes and chains are in bijection: a class-switching model with two classes
 * folded into ONE chain takes their multi-class branch against a single-chain
 * document, and reads the $U column as the residence time (the JAR's length
 * guard turns the same case into an all-zero table instead).
 */
inline bool parse_row(const std::string& line, std::size_t nchains, ParsedRow* out) {
    if (line.find(',') == std::string::npos) return false;
    if (line.find('$') != std::string::npos) return false;  // the header
    std::string s;
    for (std::size_t i = 0; i < line.size(); ++i)
        if (line[i] != ' ' && line[i] != '\t' && line[i] != '\r') s += line[i];

    std::vector<std::string> parts;
    std::size_t b = 0;
    while (true) {
        const std::size_t j = s.find(',', b);
        parts.push_back(s.substr(b, j == std::string::npos ? j : j - b));
        if (j == std::string::npos) break;
        b = j + 1;
    }
    const std::size_t stride = nchains == 1 ? 1 : nchains + 1;
    if (parts.size() < 1 + 4 * stride) return false;

    out->station = parts[0];
    out->Q.assign(nchains, 0.0);
    out->W.assign(nchains, 0.0);
    out->U.assign(nchains, 0.0);
    out->Tp.assign(nchains, 0.0);
    std::size_t ptr = 1;
    for (std::size_t c = 0; c < nchains; ++c) out->Q[c] = parse_double_or_zero(parts[ptr + c]);
    ptr += stride;
    for (std::size_t c = 0; c < nchains; ++c) out->W[c] = parse_double_or_zero(parts[ptr + c]);
    ptr += stride;
    for (std::size_t c = 0; c < nchains; ++c) out->U[c] = parse_double_or_zero(parts[ptr + c]);
    ptr += stride;
    for (std::size_t c = 0; c < nchains; ++c) out->Tp[c] = parse_double_or_zero(parts[ptr + c]);
    return true;
}

}  // namespace detail

/**
 * The gate `runAnalyzer.m` reaches through `runAnalyzerChecks`, narrowed to what
 * the JMVA document can actually carry.
 *
 * The reference's feature set admits nodes -- a Cache, a Place, a Transition --
 * that `writeJMVA` then drops on the floor, and the tool answers a SMALLER model
 * than the caller handed it with no indication that it did. A station the
 * document cannot express is refused by name here instead.
 */
template <class T>
void check_supported(const qn::NetworkStruct<T>& L) {
    for (std::size_t i = 0; i < L.nstations; ++i) {
        const qn::NodeType nt = L.stations[i].nodetype;
        if (nt == qn::NodeType::Queue || nt == qn::NodeType::Delay ||
            nt == qn::NodeType::Source)
            continue;
        throw UnsupportedError(
            "SolverQNS: station '" + L.nodes[L.station_to_node[i] - 1].name +
            "' is not a Queue, a Delay or a Source, and the JMVA document qnsolver reads "
            "carries no other station type");
    }
    if (L.has_priorities())
        throw UnsupportedError(
            "SolverQNS: class priorities are outside the product-form envelope qnsolver "
            "evaluates");
}

/**
 * Port of `@@SolverQNS/runAnalyzer.m` and `solver_qns.m`.
 *
 * @param L   the refreshed struct
 * @param opt the method, the multiserver rule and the JMVA sample cap
 */
template <class T>
mva::AvgResult<T> solver_qns_run_analyzer(const qn::NetworkStruct<T>& L, const QnsOptions& opt) {
    check_method(opt.method);
    // NOTHING under the QNS tree reads `sn.cap` or `sn.classcap`, on EITHER of
    // the two routes below: the JMVA document is read by `qnsolver`, whose
    // MVA-family algorithms have no representation of a finite buffer, and the
    // QN2LQN route hands the model to LQNS, which has none either. A capped
    // station was therefore solved as an unbounded one and the table reported
    // the unconstrained answer under this solver's name. There is no
    // feature-registry name for plain capacity, hence the structural test --
    // `SolverMVA`, `SolverNC`, `SolverAG` and `SolverFLD` gate the same way,
    // through this same helper. It sits HERE rather than in `check_supported`,
    // which is deliberately the narrower JMVA-document gate and is skipped on
    // the layered route.
    qn::check_binding_capacity("SolverQNS", L);

    // IMMEDIATE FEEDBACK keeps a self-looping job on its server instead of
    // re-queueing it, and neither path below can state that: the JMVA document
    // `qnsolver` reads carries a mean demand and a visit count per chain, and
    // the LQN `qn2lqn` writes turns the routing into OR-fork precedences of
    // pseudo-activities on the reference task, where a repeated visit is a NEW
    // CALL. Either would answer for re-queueing under this solver's name.
    // SolverMVA warns and solves on the same property; here it is a refusal
    // because the two conversions cannot represent the self-loop at all.
    // Mirrors matlab/src/solvers/wrappers/QNS/qns_immfeed_refusal.m.
    if (L.has_immediate_feedback())
        throw UnsupportedError(
            "SolverQNS does not support immediate feedback (sn.immfeed): neither the JMVA "
            "document qnsolver reads nor the LQN QN2LQN writes can keep a self-looping job on "
            "its server. Use SolverCTMC or SolverSSA, whose state space carries the self-loop.");

    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;

    bool has_open = false;
    for (std::size_t k = 0; k < K; ++k)
        if (!std::isfinite(L.classes[k].population)) has_open = true;

    if (!L.has_product_form() && !has_open) {
        if (L.has_priorities())
            throw UnsupportedError(
                "SolverQNS: class priorities are outside the QN2LQN conversion envelope");

        const lqn::LqnModel<T> layered = io::qn2lqn(L);
        lqns::LqnsOptions lo;
        lo.multiserver = resolve_multiserver(opt);
        lo.samples = static_cast<double>(opt.samples);
        lo.keep = opt.keep;
        lo.timeout_seconds = opt.timeout;
        lqns::SolverLQNS<T> solver(layered, lo);
        const lqns::LqnsSolution<T> ls = solver.get_ensemble_avg();
        const lqn::LqnStruct<T>& lsn = solver.get_struct();
        const T zero = num_traits<T>::from_int(0);
        Matrix<T> Q(M, K, zero), U(M, K, zero), R(M, K, zero), Tp(M, K, zero);

        const auto element = [&](const std::string& name, lqn::LqnElement kind) {
            for (std::size_t idx = 1; idx <= lsn.nidx; ++idx)
                if (lsn.type[idx] == kind && lsn.names[idx] == name) return idx;
            return std::size_t(0);
        };
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t node = L.station_to_node[i] - 1;
            const qn::NodeType type = L.nodes[node].nodetype;
            if (type != qn::NodeType::Queue && type != qn::NodeType::Delay) continue;
            for (std::size_t r = 0; r < K; ++r) {
                const std::string suffix = std::to_string(node + 1) + "_" +
                                           std::to_string(r + 1);
                const std::size_t a = element("Q" + suffix, lqn::LqnElement::ACTIVITY);
                if (a == 0) continue;  // the class does not visit this station
                if (ls.defined_Q[a]) Q(i, r) = ls.QN[a];
                if (ls.defined_U[a]) U(i, r) = ls.UN[a];
                if (ls.defined_R[a]) R(i, r) = ls.RN[a];
                if (ls.defined_T[a]) Tp(i, r) = ls.TN[a];

                // Some LQNS releases omit activity utilization/throughput but
                // report the bound entry phase.  The JAR carries the same
                // fallback, and it changes no value when the activity row is
                // present.
                const std::size_t e = element("E" + suffix, lqn::LqnElement::ENTRY);
                if (e != 0 && !ls.defined_U[a] && ls.defined_U[e]) U(i, r) = ls.UN[e];
                if (e != 0 && !ls.defined_T[a] && ls.defined_T[e]) Tp(i, r) = ls.TN[e];

                const double servers = L.stations[i].nservers;
                if (std::isfinite(servers) && servers > 0.0)
                    U(i, r) = T(U(i, r) / num_traits<T>::from_double(servers));
            }
        }

        mva::AvgResult<T> out;
        out.QN = mva::filter_metric(L, Q, mva::MetricKind::QLen, nullptr);
        out.UN = mva::filter_metric(L, U, mva::MetricKind::Util, nullptr);
        out.RN = mva::filter_metric(L, R, mva::MetricKind::RespT, nullptr);
        out.TN = mva::filter_metric(L, Tp, mva::MetricKind::Tput, nullptr);
        out.AN = mva::filter_metric(L, mva::sn_get_arvr_from_tput(L, out.TN),
                                    mva::MetricKind::ArvR, nullptr);
        // runAnalyzer.m ultimately passes an empty residence-time table to
        // setAvgResults, which derives it from response time and visits.
        out.WN = mva::filter_metric(L, mva::sn_get_residt_from_respt(L, out.RN),
                                    mva::MetricKind::ResidT, nullptr);
        out.CN.clear();
        out.XN.clear();
        out.method = opt.method;
        const std::string actual = resolve_multiserver(opt);
        out.actualmethod = opt.method == "default" ? ("default/" + actual) : actual;
        out.iter = ls.iterations;
        return out;
    }

    // This is the narrower capability gate of the JMVA document.  The layered
    // path above legitimately contains routing and Join nodes that JMVA cannot
    // encode, so applying it before dispatch would reject the very models
    // QN2LQN exists to convert.
    check_supported(L);

    const std::string ms = resolve_multiserver(opt);

    // The reference passes -m only when a multiserver station is present, so a
    // single-server model answers identically under every method name.
    bool has_multiserver = false;
    for (std::size_t i = 0; i < M; ++i) {
        const double c = L.stations[i].nservers;
        if (std::isfinite(c) && c > 1.0) has_multiserver = true;
    }

    // THE GATE BELONGS INSIDE THE MULTISERVER BRANCH, because that is the only
    // branch the approximation reaches. Without a multiserver station the
    // reference emits no -m at all and answers under the caller's method name,
    // so refusing 'suri' there would refuse a model the reference solves; with
    // one, `solver_qns.m` falls off its switch and leaves `cmd` unassigned,
    // which is an undefined-variable error rather than a diagnosis.
    if (has_multiserver && !is_qnsolver_multiserver(ms))
        throw UnsupportedError(
            "SolverQNS: the multiserver approximation '" + ms +
            "' is one LQNS offers and qnsolver does not: 'qnsolver -m' accepts conway, reiser, "
            "rolia and zhou only; suri and schmidt are available only on the non-product-form "
            "closed SolverLQNS branch");

    if (!is_available())
        throw UnsupportedError(
            "SolverQNS needs the external 'qnsolver' binary on the PATH. It ships with LQNS "
            "(http://www.sce.carleton.ca/rads/lqns/); LINE distributes no copy and runs none "
            "from a container image, because that licence forbids redistribution");

    util::TempDir tmp("qns");
    if (opt.keep) tmp.keep();
    const std::string model_file = tmp.file("model.jmva");
    const std::string result_file = tmp.file("result.jmva");
    line::util::LineConsole::step("writing the JMVA model file");
    io::write_jmva(L, model_file, opt.method, opt.samples);

    std::vector<std::string> argv = detail::qnsolver_argv();
    // `-l` IS `--linearizer`, AND IT TAKES NO ARGUMENT: `qnsolver --help` lists
    // it beside `-e, --exact-mva` and `-s, --schweitzer` as a solver selector,
    // and the input file is POSITIONAL. So this line does not mean "load the
    // model" -- it selects Linearizer and lets the model file fall through as
    // the positional argument. All four codebases emit it, so they agree, and
    // that agreement is on LINEARIZER results rather than on the exact MVA
    // qnsolver runs by default. Left as it is deliberately: changing it changes
    // every QNS number in every codebase at once, which is the user's call.
    argv.push_back("-l");
    argv.push_back(model_file);
    if (has_multiserver) argv.push_back("-m" + ms);
    argv.push_back("-o");
    argv.push_back(result_file);

    line::util::LineConsole::step("running the qnsolver binary as a subprocess");
    const util::ProcResult pr = util::capture(argv, opt.timeout);
    if (pr.timedOut)
        throw NumericError("SolverQNS: qnsolver did not finish within " +
                           std::to_string(opt.timeout) + "s and was killed");
    if (pr.exitCode != 0)
        throw NumericError("SolverQNS: qnsolver exited with code " +
                           std::to_string(pr.exitCode) +
                           (pr.out.empty() ? std::string() : ("\n" + pr.out)));

    // ---- parse the chain-level table -------------------------------------
    line::util::LineConsole::step("parsing the qnsolver output");
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Qchain(M, C, zero), Uchain(M, C, zero), Wchain(M, C, zero), Tchain(M, C, zero);
    std::ifstream in(result_file.c_str());
    if (!in)
        // qnsolver exits 0 on a parse error, so the exit-code branch above never fires and its
        // own message is the only evidence of what it rejected. Carry it here or a build that
        // cannot read an <ldstation> reads as an unexplained missing file.
        throw NumericError("SolverQNS: qnsolver wrote no result file at '" + result_file + "'" +
                           (pr.out.empty() ? std::string() : ("\nqnsolver said: " + pr.out)));
    std::string line;
    std::size_t nrows = 0;
    while (std::getline(in, line)) {
        detail::ParsedRow row;
        if (!detail::parse_row(line, C, &row)) continue;
        std::size_t idx = M;
        for (std::size_t i = 0; i < M; ++i)
            if (L.nodes[L.station_to_node[i] - 1].name == row.station) {
                idx = i;
                break;
            }
        if (idx == M) continue;  // a station qnsolver names and the model does not
        for (std::size_t c = 0; c < C; ++c) {
            Qchain(idx, c) = num_traits<T>::from_double(row.Q[c]);
            Wchain(idx, c) = num_traits<T>::from_double(row.W[c]);
            Uchain(idx, c) = num_traits<T>::from_double(row.U[c]);
            Tchain(idx, c) = num_traits<T>::from_double(row.Tp[c]);
        }
        ++nrows;
    }
    if (nrows == 0)
        throw NumericError(
            "SolverQNS: qnsolver produced no station rows this model recognises; its output "
            "names no station of the model");

    const mva::ChainDemands<T> d = mva::sn_get_demands_chain(L);

    // ---- chain throughput at the reference station -----------------------
    // THE REFERENCE STATION IS READ PER CHAIN, through the chain's first class.
    // MATLAB and the JAR index `sn.refstat` -- a per-CLASS array -- with the
    // chain number, which lands on an unrelated class's reference station as
    // soon as class switching makes the two index spaces differ.
    std::vector<T> Xchain(C, zero);
    for (std::size_t c = 0; c < C; ++c) {
        const std::size_t rstat = L.classes[L.inchain[c][0] - 1].refstat;
        if (Tchain(rstat - 1, c) > zero) {
            Xchain[c] = Tchain(rstat - 1, c);
            continue;
        }
        // An open chain's reference station is the Source, which the JMVA
        // document does not carry, so its row is zero; recover X from any
        // station whose visit count is known.
        for (std::size_t i = 0; i < M; ++i)
            if (d.Vchain(i, c) > zero && Tchain(i, c) > zero) {
                Xchain[c] = T(Tchain(i, c) / d.Vchain(i, c));
                break;
            }
    }

    // `Rchain = Wchain` of solver_qns.m: qnsolver's $R column is already the
    // per-visit residence at the station, so it is not divided by the visits.
    Matrix<T> Rchain = Wchain;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c)
            if (std::isnan(num_traits<T>::to_double(Rchain(i, c)))) Rchain(i, c) = zero;

    // Utilization comes back summed over the servers of a multiserver station,
    // which the ld encoding of writeJMVA turns it into; LINE reports it per
    // server.
    for (std::size_t i = 0; i < M; ++i) {
        const double c = L.stations[i].nservers;
        if (!std::isfinite(c)) continue;
        for (std::size_t j = 0; j < C; ++j)
            Uchain(i, j) = T(Uchain(i, j) / num_traits<T>::from_double(c));
    }

    const mva::ClassResults<T> cr =
        mva::sn_deaggregate_chain_results(L, d, Qchain, Uchain, Rchain, Tchain, Xchain);

    mva::AvgResult<T> out;
    out.QN = mva::filter_metric(L, cr.Q, mva::MetricKind::QLen, nullptr);
    out.UN = mva::filter_metric(L, cr.U, mva::MetricKind::Util, nullptr);
    out.RN = mva::filter_metric(L, cr.R, mva::MetricKind::RespT, nullptr);
    out.TN = mva::filter_metric(L, cr.Tp, mva::MetricKind::Tput, nullptr);
    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source)
            for (std::size_t k = 0; k < K; ++k) srcmask[i][k] = true;
    out.AN = mva::filter_metric(L, mva::sn_get_arvr_from_tput(L, out.TN), mva::MetricKind::ArvR,
                                &srcmask);
    // MATLAB passes [] for the residence time and lets `getAvg` derive it from
    // the response time and the visits, which is what this helper is; the JAR
    // instead sets WN = RN, and the two agree only where the visit count is one.
    out.WN = mva::filter_metric(L, mva::sn_get_residt_from_respt(L, out.RN),
                                mva::MetricKind::ResidT, nullptr);
    out.CN = cr.C;
    out.XN = cr.X;
    out.method = opt.method;
    // `default` is reported as 'default/<what ran>', the convention of
    // runAnalyzer.m. Without a multiserver station no -m flag is passed, and the
    // reference then leaves the label at the tool's own default.
    const std::string ran = has_multiserver ? ms : std::string("default");
    out.actualmethod = opt.method == "default" ? ("default/" + ran) : ran;
    out.iter = 0;
    return out;
}

}  // namespace qns
}  // namespace line

#endif  // LINE_SOLVERS_WRAPPERS_QNS_SOLVER_QNS_H
