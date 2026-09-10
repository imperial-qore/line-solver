/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_WRAPPERS_LQNS_SOLVER_LQNS_H
#define LINE_SOLVERS_WRAPPERS_LQNS_SOLVER_LQNS_H

/**
 * SolverLQNS: the layered model solved by the external lqns / lqsim binaries.
 *
 * Port of matlab/src/solvers/wrappers/LQNS/@@SolverLQNS (SolverLQNS.m,
 * runAnalyzer.m, parseXMLResults.m, getEnsembleAvg.m, runRemoteLQNS.m), of
 * jline.solvers.wrappers.lqns.SolverLQNS and of the native-Python twin. It is a
 * WRAPPER, not an engine: nothing here computes a queueing result. The model is
 * written as .lqnx (lqn_writer.h), a binary is run over it, and the .lqxo it
 * leaves behind is read back.
 *
 * LINE SHIPS NO LQNS BINARY. Its licence is an evaluation agreement that
 * forbids redistribution, so the binary is one the user installed; a missing
 * one is reported by name, with the two ways to obtain an answer anyway (an
 * install, or a host already running the REST service), rather than as a failed
 * exec.
 *
 * WHAT THE COLUMNS MEAN, and why they are not the ones lqns prints. The
 * reference's getEnsembleAvg permutes them so that a layered result reads the
 * same whichever solver produced it:
 *
 *   QLen  <- the element's UTILIZATION      (a task's utilization IS its
 *                                            mean number in service)
 *   Util  <- its PROCESSOR utilization: lqns and SolverLN both sum it over
 *            the host's servers, so no rescaling applies. Verbatim for hosts,
 *            tasks and activities; for an entry it is aggregated over the
 *            activity graph, which lqns itself reports as 0 in the
 *            activity-graph form
 *   RespT <- its PHASE 1 SERVICE TIME
 *   Tput  <- its throughput
 *   ResidT, ArvR: lqns computes neither, so they stay undefined rather than
 *            being filled with a zero that would read as a computed value
 *
 * THE .lqxo IS THE ONLY OUTPUT READ. lqns also prints a human-readable report;
 * parsing that instead would tie the port to a print format that changes
 * between releases, and it carries fewer digits than the XML.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "json.hpp"
#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/lqn/lqn_writer.h"
#include "line/solvers/wrappers/lqns/lqns_probe.h"
#include "line/util/error.h"
#include "line/util/http.h"
#include "line/util/subprocess.h"
#include "line/util/tempdir.h"
#include "line/util/xml.h"

namespace line {
namespace lqns {

/** Knobs of the wrapper, the subset of SolverOptions that reaches lqns. */
struct LqnsOptions {
    /**
     * default | lqns | srvn | exactmva | srvn.exactmva | sim | lqsim | lqnsdefault
     *
     * `sim` and `lqsim` run the SIMULATOR and are the only stochastic ones;
     * `lqnsdefault` is lqns with no pragma at all, which is a different fixed
     * point from `default` and not a synonym for it.
     */
    std::string method = "default";
    /** conway | rolia | zhou | suri | reiser | schmidt | default (= rolia). */
    std::string multiserver = "default";
    /** lqsim run length, `-A`; not positive leaves lqsim's own default. */
    double samples = 10000.0;
    bool verbose = false;
    /** Keep the working directory, model and result file after the run. */
    bool keep = false;
    /** Deadline for the child, in seconds; not positive waits indefinitely. */
    int timeout_seconds = 0;
    /** Solve on a host running lqns-rest instead of on this machine. */
    bool remote = false;
    std::string remote_url = "http://localhost:8080";
};

/**
 * Everything the .lqxo carries, indexed as the struct is.
 *
 * The element vectors are (nidx+1) with slot 0 unused and `callwaiting` is
 * (ncalls+1); an entry lqns did not report stays NaN, which is how the
 * reference leaves it and is what tells a caller "not reported" apart from
 * "reported as zero".
 */
struct LqnsRawAvg {
    std::vector<double> util, phase1util, phase2util;
    std::vector<double> phase1svct, phase2svct;
    std::vector<double> tput, procwaiting, procutil;
    std::vector<double> callwaiting;
    int iterations = 0;
};

/** The six measures, on the element index space, with a defined mask each. */
template <class T>
struct LqnsSolution {
    std::vector<T> QN, UN, RN, TN, AN, WN;
    std::vector<bool> defined_Q, defined_U, defined_R, defined_T, defined_A, defined_W;
    int iterations = 0;
};

namespace detail {

/** NaN as the .lqxo reader uses it: "lqns did not report this". */
inline double nan_value() { return std::numeric_limits<double>::quiet_NaN(); }

/** An attribute as a double, NaN when absent or unparsable (str2double). */
inline double attr_num(const xml::Element* e, const char* key) {
    const std::string s = e->attr(key);
    if (s.empty()) return nan_value();
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str()) return nan_value();
    return v;
}

/**
 * Snap a value within CoarseTol of an exact tenth onto it.
 *
 * lqns reports a quantity that is analytically a multiple of 0.1 with a few
 * digits of iteration noise, and the reference removes that noise in
 * getAvgTable before tabulating. It belongs to the TABLE and not to getAvg: the
 * raw fixed point is what a parity comparison should see.
 */
inline double snap_to_tenth(double v) {
    if (!std::isfinite(v)) return v;
    const double scaled = v * 10.0;
    const double snapped = std::floor(scaled + 0.5);
    if (std::fabs(scaled - snapped) < lang::GlobalConstants::CoarseTol * scaled)
        return snapped / 10.0;
    return v;
}

/**
 * A private working directory, the C++ lineTempName.
 *
 * Delegated rather than rolled again here: LINE_WORKSPACE_ROOT has to relocate
 * the staged .lqnx, since a containerized lqns bind-mounts that root alone and
 * sees nothing left under TMPDIR.
 */
inline std::string make_temp_dir(const std::string& tag) { return util::make_temp_dir(tag); }

/** Remove a directory this wrapper created, with every file it holds. */
inline void remove_temp_dir(const std::string& dir) {
    DIR* d = ::opendir(dir.c_str());
    if (d != nullptr) {
        for (struct dirent* ent = ::readdir(d); ent != nullptr; ent = ::readdir(d)) {
            const std::string name(ent->d_name);
            if (name == "." || name == "..") continue;
            ::unlink((dir + "/" + name).c_str());
        }
        ::closedir(d);
    }
    ::rmdir(dir.c_str());
}

}  // namespace detail

/**
 * The layered model solved by lqns or lqsim.
 *
 * Constructed over the INTERMEDIATE model (lqn_reader.h), not over the struct,
 * because the document handed to the binary has to carry the precedence blocks
 * and reply entries that getStruct flattens away.
 */
template <class T>
class SolverLQNS {
public:
    SolverLQNS(const lqn::LqnModel<T>& model, const LqnsOptions& options = LqnsOptions())
        : model_(model), opt_(options), sn_(lqn::lqn_finalize(model)) {
        if (!is_available() && !opt_.remote)
            throw UnsupportedError(
                "SolverLQNS requires the lqns and lqsim commands on the system path.\n"
                "Obtain them from their authors at: http://www.sce.carleton.ca/rads/lqns/\n"
                "LINE ships no LQNS binary and does not redistribute one.\n\n"
                "Alternatively, point LINE at a host that already runs LQNS by setting\n"
                "LqnsOptions::remote and LqnsOptions::remote_url.");
        const std::string m = opt_.method.empty() ? std::string("default") : opt_.method;
        const std::vector<std::string> valid = list_valid_methods();
        if (std::find(valid.begin(), valid.end(), m) == valid.end())
            throw InputError("SolverLQNS: '" + m +
                             "' is not a method of this solver; it takes default, lqns, srvn, "
                             "exactmva, srvn.exactmva, sim, lqsim and lqnsdefault");
        opt_.method = m;
        multiserver_pragma(opt_.multiserver);  // refuses an unknown policy here, not mid-run
    }

    /** True when the native lqns command answers on this machine, any release. */
    static bool has_local_binary() { return !lqns_version().empty(); }

    /** True when a local binary is installed AND is 6.0 or greater. */
    static bool is_available() { return lqns_is_available(); }

    /** The version banner of the local binary, empty when there is none. */
    static std::string version() { return lqns_version(); }

    static std::vector<std::string> list_valid_methods() {
        std::vector<std::string> m;
        m.push_back("default");
        m.push_back("lqns");
        m.push_back("srvn");
        m.push_back("exactmva");
        m.push_back("srvn.exactmva");
        m.push_back("sim");
        m.push_back("lqsim");
        m.push_back("lqnsdefault");
        return m;
    }

    /** Only the lqsim methods draw random numbers. */
    static bool is_stochastic_method(const std::string& method) {
        return method == "sim" || method == "lqsim";
    }

    /**
     * The per-layer feature set, as SolverLQNS.supports declares it.
     *
     * A layered model reaches lqns as a whole, so this is the set each LAYER
     * must fall inside: the product-form station kinds and the four service
     * distributions the LQN schema can carry.
     */
    static std::vector<std::string> feature_set() {
        const char* names[] = {"Sink",          "Source",        "Queue",
                               "Coxian",        "Erlang",        "Exp",
                               "HyperExp",      "Buffer",        "Server",
                               "JobSink",       "RandomSource",  "ServiceTunnel",
                               "SchedStrategy_PS", "SchedStrategy_FCFS", "ClosedClass"};
        return std::vector<std::string>(names, names + sizeof(names) / sizeof(names[0]));
    }

    const lqn::LqnStruct<T>& get_struct() const { return sn_; }
    const LqnsOptions& options() const { return opt_; }
    /** Wall-clock seconds of the last run, the reference's `runtime`. */
    double runtime() const { return runtime_; }
    int iterations() const { return raw_.iterations; }
    /** Everything the .lqxo carried, before the column permutation. */
    const LqnsRawAvg& raw_avg() {
        run_analyzer_once();
        return raw_;
    }

    /**
     * Run the binary and read its result back.
     *
     * Re-running is what the reference does on every getEnsembleAvg call, but
     * repeating a lqsim run would silently change the answer between two reads
     * of the same solver object, so the result is computed once and kept.
     */
    void run_analyzer() {
        const std::clock_t t0 = std::clock();
        const std::string dir = detail::make_temp_dir("lqns");
        const std::string stem = dir + "/model";
        const std::string modelfile = stem + ".lqnx";
        const std::string resultfile = stem + ".lqxo";
        try {
            line::util::LineConsole::step("writing the LQN model to %s", modelfile.c_str());
            const lqn::LqnWriteReport rep = lqn::write_lqnx(model_, modelfile, "LQN");
            if (opt_.verbose)
                for (std::size_t i = 0; i < rep.dropped.size(); ++i)
                    std::fprintf(stderr, "SolverLQNS: %s\n", rep.dropped[i].c_str());

            if (opt_.remote) run_remote(modelfile, resultfile);
            else run_local(modelfile);

            line::util::LineConsole::step("parsing the lqns XML results");
            parse_lqxo(resultfile);
        } catch (...) {
            if (!opt_.keep) detail::remove_temp_dir(dir);
            throw;
        }
        if (!opt_.keep) detail::remove_temp_dir(dir);
        else if (opt_.verbose) std::fprintf(stderr, "SolverLQNS: files kept in %s\n", dir.c_str());
        runtime_ = double(std::clock() - t0) / double(CLOCKS_PER_SEC);
        solved_ = true;
    }

    /**
     * The six measures on the element index space.
     *
     * Port of getEnsembleAvg.m, including the host-multiplicity rescaling of
     * the utilization column.
     */
    LqnsSolution<T> get_ensemble_avg() {
        run_analyzer_once();
        const std::size_t n = sn_.nidx;
        LqnsSolution<T> s;
        s.iterations = raw_.iterations;
        auto fill = [&](std::vector<T>& v, std::vector<bool>& d, const std::vector<double>& src) {
            v.assign(n + 1, num_traits<T>::from_int(0));
            d.assign(n + 1, false);
            for (std::size_t i = 1; i <= n; ++i) {
                if (std::isnan(src[i])) continue;
                v[i] = num_traits<T>::from_double(src[i]);
                d[i] = true;
            }
        };
        // Both lqns and SolverLN report the processor utilization summed over
        // the host's servers, so no rescaling applies; the entry rows were
        // aggregated over the activity graph by aggregate_entry_procutil.
        fill(s.QN, s.defined_Q, raw_.util);
        fill(s.UN, s.defined_U, raw_.procutil);
        fill(s.RN, s.defined_R, raw_.phase1svct);
        fill(s.TN, s.defined_T, raw_.tput);
        // lqns computes neither a residence time nor an arrival rate per
        // element; leaving them undefined is the whole of the claim.
        s.AN.assign(n + 1, num_traits<T>::from_int(0));
        s.WN.assign(n + 1, num_traits<T>::from_int(0));
        s.defined_A.assign(n + 1, false);
        s.defined_W.assign(n + 1, false);
        return s;
    }

private:
    void run_analyzer_once() {
        if (!solved_) run_analyzer();
    }

    /** The `-Pmultiserver=` argument of a policy name, empty when there is none. */
    static std::string multiserver_pragma(const std::string& policy) {
        if (policy.empty() || policy == "none") return std::string();
        if (policy == "default") return "-Pmultiserver=rolia";
        if (policy == "conway" || policy == "rolia" || policy == "zhou" || policy == "suri" ||
            policy == "reiser" || policy == "schmidt")
            return "-Pmultiserver=" + policy;
        throw InputError("SolverLQNS: '" + policy +
                         "' is not a multiserver policy of lqns; it takes conway, rolia, zhou, "
                         "suri, reiser, schmidt and default");
    }

    /**
     * The command line, as an argv vector.
     *
     * NO SHELL SEES THESE ARGUMENTS. The reference builds one string and hands
     * it to system(); a model path holding a space would be word-split there,
     * and this port's working directory comes from TMPDIR, which the caller
     * controls.
     *
     * `env -u LD_LIBRARY_PATH`, which the MATLAB wrapper prefixes, is NOT
     * reproduced: it exists because MATLAB injects its own libstdc++ ahead of
     * the system one and lqns then fails to find a GLIBCXX symbol. A plain C++
     * process has no such injection, and stripping the variable here would
     * instead break a user who set it to reach their own lqns.
     */
    std::vector<std::string> build_argv(const std::string& modelfile) const {
        const bool sim = is_stochastic_method(opt_.method);
        std::vector<std::string> a;
        a.push_back(sim ? "lqsim" : "lqns");
        if (!opt_.verbose) {
            a.push_back("-a");  // no advisories
            a.push_back("-w");  // no warnings
        }
        // The simulator has no MVA to configure, so the multiserver pragma is
        // not passed to it, exactly as the reference declines to.
        if (!sim) {
            const std::string ms = multiserver_pragma(opt_.multiserver);
            if (!ms.empty()) a.push_back(ms);
        }
        if (opt_.method == "srvn" || opt_.method == "srvn.exactmva") a.push_back("-Playering=srvn");
        if (opt_.method == "exactmva" || opt_.method == "srvn.exactmva") a.push_back("-Pmva=exact");
        if (sim && opt_.samples > 0) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.0f", opt_.samples);
            a.push_back("-A");
            a.push_back(buf);
        }
        // `lqnsdefault` is lqns with NO pragma: a model that loses messages
        // stops the run there, which is the binary's own default and a
        // different answer from the one the other methods ask for.
        if (opt_.method != "lqnsdefault") a.push_back("-Pstop-on-message-loss=false");
        a.push_back("-x");  // XML result, the .lqxo this wrapper reads
        a.push_back(modelfile);
        return a;
    }

    void run_local(const std::string& modelfile) const {
        const std::vector<std::string> argv = build_argv(modelfile);
        if (opt_.verbose) {
            std::string cmd;
            for (std::size_t i = 0; i < argv.size(); ++i) cmd += (i ? " " : "") + argv[i];
            std::fprintf(stderr, "SolverLQNS command: %s\n", cmd.c_str());
        }
        line::util::LineConsole::step("running the lqns binary as a subprocess");
        const util::ProcResult r = util::capture(argv, opt_.timeout_seconds);
        if (r.timedOut)
            throw Error("SolverLQNS: " + argv[0] + " did not finish within " +
                        std::to_string(opt_.timeout_seconds) + "s and was killed");
        if (r.exitCode < 0)
            throw UnsupportedError("SolverLQNS: could not run '" + argv[0] +
                                   "'; LINE ships no LQNS binary, so it must be installed and on "
                                   "the system path");
        if (opt_.verbose && !r.out.empty()) std::fprintf(stderr, "%s", r.out.c_str());
    }

    /**
     * The lqns-rest protocol: one POST carrying the whole document.
     *
     * Twin of runRemoteLQNS.m and of the native-Python `_run_remote_lqns`. The
     * service answers with the .lqxo text, which is written where the local run
     * would have left it so that ONE parser serves both paths.
     */
    void run_remote(const std::string& modelfile, const std::string& resultfile) const {
        std::string content;
        {
            std::ifstream in(modelfile.c_str());
            if (!in) throw Error("SolverLQNS: cannot read the model it just wrote: " + modelfile);
            content.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
        }
        const bool sim = is_stochastic_method(opt_.method);
        std::string base = opt_.remote_url;
        while (!base.empty() && base[base.size() - 1] == '/') base.erase(base.size() - 1);
        const std::string url = base + (sim ? "/api/v1/solve/lqsim" : "/api/v1/solve/lqns");

        nlohmann::json req;
        req["model"]["content"] = content;
        req["model"]["base64"] = false;
        req["options"]["include_raw_output"] = true;
        if (sim) {
            req["options"]["blocks"] = 30;
            if (opt_.samples > 0) req["options"]["run_time"] = opt_.samples;
        } else {
            std::string ms = opt_.multiserver.empty() || opt_.multiserver == "default"
                                 ? std::string("rolia")
                                 : opt_.multiserver;
            req["options"]["pragmas"]["multiserver"] = ms;
            req["options"]["pragmas"]["stop_on_message_loss"] = false;
            if (opt_.method == "srvn" || opt_.method == "srvn.exactmva")
                req["options"]["pragmas"]["layering"] = "srvn";
            if (opt_.method == "exactmva" || opt_.method == "srvn.exactmva")
                req["options"]["pragmas"]["mva"] = "exact";
        }

        const int timeout_ms =
            (opt_.timeout_seconds > 0 ? opt_.timeout_seconds : 300) * 1000;
        const http::Response resp = http::post_json(url, req.dump(), timeout_ms);
        if (resp.status < 200 || resp.status >= 300)
            throw Error("SolverLQNS: remote LQNS at " + url + " returned HTTP " +
                        std::to_string(resp.status) + ": " + resp.body);
        nlohmann::json out;
        try {
            out = nlohmann::json::parse(resp.body);
        } catch (const std::exception& e) {
            throw Error(std::string("SolverLQNS: remote LQNS returned a body that is not JSON: ") +
                        e.what());
        }
        const std::string status = out.value("status", std::string());
        if (status == "error" || status == "failed")
            throw Error("SolverLQNS: remote solver returned an error: " +
                        out.value("error", std::string("unknown error")));
        std::string lqxo;
        if (out.contains("raw_output") && out["raw_output"].is_object())
            lqxo = out["raw_output"].value("lqxo", std::string());
        if (lqxo.empty())
            throw Error("SolverLQNS: remote solver did not return an LQXO document");
        std::ofstream f(resultfile.c_str());
        if (!f) throw Error("SolverLQNS: cannot write the remote result to " + resultfile);
        f << lqxo;
    }

    /**
     * Element index of a name OF A GIVEN KIND, 0 when we hold no such element.
     *
     * THE KIND IS PART OF THE KEY, and has to be. A LINE-generated layered model
     * routinely gives a processor, its task and that task's entry the SAME name
     * (`c0` in the randomLQN corpus, and throughout lqn_ofbiz), and `lqn.names`
     * holds all three. The reference looks the name up in that flat list with
     * `findstring`, which returns EVERY match, and then assigns the result to
     * all of them -- so a processor row ends up carrying its entry's numbers,
     * whichever element was written last. The .lqxo says which kind it is
     * describing, in the tag being read, so the ambiguity does not have to exist
     * here and this port does not reproduce it. On a model whose names are
     * unique the two agree element for element.
     */
    std::size_t index_of(const std::string& name, lang::LqnElement kind) const {
        const std::map<std::pair<int, std::string>, std::size_t>::const_iterator it =
            byname_.find(std::make_pair(static_cast<int>(kind), name));
        return it == byname_.end() ? 0 : it->second;
    }

    std::size_t call_index_of(const std::string& name) const {
        const std::map<std::string, std::size_t>::const_iterator it = bycallname_.find(name);
        return it == bycallname_.end() ? 0 : it->second;
    }

    void init_name_index() {
        if (!byname_.empty()) return;
        // First declaration wins within a kind, as `find(strcmp(...))(1)` would.
        for (std::size_t i = 1; i <= sn_.nidx; ++i)
            byname_.insert(std::make_pair(
                std::make_pair(static_cast<int>(sn_.type[i]), sn_.names[i]), i));
        for (std::size_t c = 1; c <= sn_.ncalls; ++c)
            bycallname_.insert(std::make_pair(sn_.callnames[c], c));
    }

    /**
     * Read the .lqxo, a port of parseXMLResults.m.
     *
     * The walk is processor -> task -> entry -> phase activities, then the task
     * activity graph. It is keyed on NAMES throughout, so an activity ordering
     * that differs between the document and the struct cannot mis-assign a row.
     */
    void parse_lqxo(const std::string& resultfile) {
        init_name_index();
        const std::size_t n = sn_.nidx;
        const double nan = detail::nan_value();
        raw_ = LqnsRawAvg();
        raw_.util.assign(n + 1, nan);
        raw_.phase1util.assign(n + 1, nan);
        raw_.phase2util.assign(n + 1, nan);
        raw_.phase1svct.assign(n + 1, nan);
        raw_.phase2svct.assign(n + 1, nan);
        raw_.tput.assign(n + 1, nan);
        raw_.procwaiting.assign(n + 1, nan);
        raw_.procutil.assign(n + 1, nan);
        raw_.callwaiting.assign(sn_.ncalls + 1, nan);

        std::unique_ptr<xml::Element> doc;
        try {
            doc = xml::parse_file(resultfile);
        } catch (const Error&) {
            throw Error(
                "SolverLQNS: no readable result at " + resultfile +
                "; the binary ran but wrote no .lqxo, which is what it does when it rejects the "
                "model (run with verbose to see its diagnostics)");
        }

        for (const xml::Element* sp : doc->by_tag("solver-params")) {
            const std::vector<const xml::Element*> g = sp->by_tag("result-general");
            if (!g.empty()) raw_.iterations = static_cast<int>(detail::attr_num(g[0], "iterations"));
        }

        for (const xml::Element* pe : doc->by_tag("processor")) {
            const std::size_t pidx = index_of(pe->attr("name"), lang::LqnElement::HOST);
            const std::vector<const xml::Element*> pr = pe->by_tag("result-processor");
            if (pidx && !pr.empty()) raw_.procutil[pidx] = detail::attr_num(pr[0], "utilization");

            for (const xml::Element* te : pe->by_tag("task")) {
                const std::size_t tidx = index_of(te->attr("name"), lang::LqnElement::TASK);
                const std::vector<const xml::Element*> tr = te->by_tag("result-task");
                double task_tput = nan;
                if (!tr.empty()) {
                    task_tput = detail::attr_num(tr[0], "throughput");
                    if (tidx) {
                        raw_.util[tidx] = detail::attr_num(tr[0], "utilization");
                        raw_.phase1util[tidx] = detail::attr_num(tr[0], "phase1-utilization");
                        raw_.phase2util[tidx] = detail::attr_num(tr[0], "phase2-utilization");
                        raw_.tput[tidx] = task_tput;
                        raw_.procutil[tidx] = detail::attr_num(tr[0], "proc-utilization");
                    }
                }

                for (const xml::Element* ee : te->by_tag("entry")) {
                    const std::size_t eidx = index_of(ee->attr("name"), lang::LqnElement::ENTRY);
                    const std::vector<const xml::Element*> er = ee->by_tag("result-entry");
                    double entry_tput = nan;
                    if (!er.empty()) {
                        entry_tput = detail::attr_num(er[0], "throughput");
                        if (eidx) {
                            raw_.util[eidx] = detail::attr_num(er[0], "utilization");
                            raw_.phase1util[eidx] = detail::attr_num(er[0], "phase1-utilization");
                            raw_.phase2util[eidx] = detail::attr_num(er[0], "phase2-utilization");
                            raw_.phase1svct[eidx] = detail::attr_num(er[0], "phase1-service-time");
                            raw_.phase2svct[eidx] = detail::attr_num(er[0], "phase2-service-time");
                            raw_.tput[eidx] = entry_tput;
                            raw_.procutil[eidx] = detail::attr_num(er[0], "proc-utilization");
                        }
                    }

                    // PH1PH2 form: the phase activities carry their own results
                    const std::vector<const xml::Element*> epa =
                        ee->by_tag("entry-phase-activities");
                    if (epa.empty()) continue;
                    for (const xml::Element* ae : epa[0]->by_tag("activity")) {
                        read_activity(ae, entry_tput, true);
                    }
                }

                const std::vector<const xml::Element*> tal = te->by_tag("task-activities");
                if (tal.empty()) continue;
                for (const xml::Element* ae : tal[0]->by_tag("activity")) {
                    // A synch-call of a phase activity is nested inside the same
                    // subtree in some documents; only direct children of
                    // task-activities are the graph's own activities.
                    if (ae->parent != tal[0]) continue;
                    read_activity(ae, task_tput, false);
                }
            }
        }

        aggregate_entry_procutil();
        aggregate_entry_svct();
    }

    /**
     * Phase-1 service time of an entry lqns never invoked.
     *
     * lqns omits `phase1-service-time` from `result-entry` exactly when the
     * entry's throughput is zero: nothing was served, so there is no
     * per-invocation mean to report. LINE then carried a NaN where the table
     * says an entry HAS a response time and every other solver reports one,
     * breaking the NaN mask -- see `_kb/06-solver-catalog.md`. The value is
     * taken from the activity rows, and ONLY where they are unanimous: if every
     * activity reachable from the entry reports a zero service time then every
     * aggregation law agrees on zero -- the serial sum, the branch-weighted mean
     * of an OrFork, the order statistic of an AndFork -- so the derivation does
     * not depend on which one applies.
     *
     * It is deliberately NOT generalised the way `aggregate_entry_procutil` is.
     * Utilizations add over an activity graph; response times do not. Measured
     * over the example corpus, the sum over actsof reproduces
     * `phase1-service-time` on serial chains only and misses it wherever the
     * graph branches (`lqn_workflows` `Entry`: 12.5667 reported against 8.5667
     * summed, `lqn_fork_open_arrival` `SE`: 0.841667 against 1.0), so a summed
     * fallback would answer with a number lqns contradicts. An entry whose
     * activities are unreported, absent, or not all zero keeps NaN.
     */
    void aggregate_entry_svct() {
        for (std::size_t e = 1; e <= sn_.nentries; ++e) {
            const std::size_t eidx = sn_.eshift + e;
            if (eidx >= sn_.actsof.size() || eidx >= raw_.phase1svct.size()) continue;
            if (!std::isnan(raw_.phase1svct[eidx])) continue;
            const std::vector<std::size_t>& acts = sn_.actsof[eidx];
            if (acts.empty()) continue;
            bool all_zero = true;
            for (std::size_t aidx : acts) {
                if (aidx >= raw_.phase1svct.size() || std::isnan(raw_.phase1svct[aidx]) ||
                    raw_.phase1svct[aidx] != 0.0) {
                    all_zero = false;
                    break;
                }
            }
            if (all_zero) raw_.phase1svct[eidx] = 0.0;
        }
    }

    /**
     * Processor utilization of an entry, aggregated from its activity graph.
     *
     * lqns credits host work to whichever level carries the host demand. In the
     * activity-graph form an entry declares none, so lqns reports its
     * result-entry proc-utilization as a literal 0 and the work sits on the
     * result-activity rows; the entry's value is then the sum over the
     * activities reachable from the entry within its own task, which is what
     * actsof holds. In PH1PH2 form the same sum runs over the phase activities
     * and reproduces the value lqns reports there, so no form test is needed. An
     * entry with no activities, or any activity lqns left unreported, keeps the
     * raw attribute rather than a partial sum.
     */
    void aggregate_entry_procutil() {
        for (std::size_t e = 1; e <= sn_.nentries; ++e) {
            const std::size_t eidx = sn_.eshift + e;
            if (eidx >= sn_.actsof.size()) continue;
            const std::vector<std::size_t>& acts = sn_.actsof[eidx];
            if (acts.empty()) continue;
            double sum = 0.0;
            bool complete = true;
            for (std::size_t aidx : acts) {
                if (aidx >= raw_.procutil.size() || std::isnan(raw_.procutil[aidx])) {
                    complete = false;
                    break;
                }
                sum += raw_.procutil[aidx];
            }
            if (complete) raw_.procutil[eidx] = sum;
        }
    }

    /**
     * One `<activity>` and the calls it issues.
     *
     * @param ae          the element
     * @param owner_tput  throughput of the entry (phase form) or task, used
     *                    where lqns omits the activity's own
     * @param phase_form  true inside entry-phase-activities, where lqns omits
     *                    both throughput and proc-utilization
     */
    void read_activity(const xml::Element* ae, double owner_tput, bool phase_form) {
        const std::string aname = ae->attr("name");
        const std::size_t aidx = index_of(aname, lang::LqnElement::ACTIVITY);
        const std::vector<const xml::Element*> ar = ae->by_tag("result-activity");
        if (aidx && !ar.empty()) {
            const xml::Element* r = ar[0];
            raw_.util[aidx] = detail::attr_num(r, "utilization");
            raw_.phase1svct[aidx] = detail::attr_num(r, "service-time");
            raw_.procwaiting[aidx] = detail::attr_num(r, "proc-waiting");
            const double t = detail::attr_num(r, "throughput");
            // Each phase executes once per entry invocation, so an omitted
            // throughput there IS the entry's; filling it with zero would say
            // the activity never runs.
            raw_.tput[aidx] = std::isnan(t) && phase_form ? owner_tput : t;
            const double pu = detail::attr_num(r, "proc-utilization");
            if (!std::isnan(pu)) {
                raw_.procutil[aidx] = pu;
            } else if (phase_form) {
                const double hd = detail::attr_num(ae, "host-demand-mean");
                if (!std::isnan(owner_tput) && !std::isnan(hd)) raw_.procutil[aidx] = owner_tput * hd;
            }
        }
        if (!aidx) return;
        read_calls(ae, aname, "synch-call", "=>");
        read_calls(ae, aname, "asynch-call", "->");
    }

    void read_calls(const xml::Element* ae, const std::string& aname, const char* tag,
                    const char* arrow) {
        for (const xml::Element* ce : ae->by_tag(tag)) {
            const std::size_t cidx = call_index_of(aname + arrow + ce->attr("dest"));
            const std::vector<const xml::Element*> cr = ce->by_tag("result-call");
            if (cidx && !cr.empty()) raw_.callwaiting[cidx] = detail::attr_num(cr[0], "waiting");
        }
    }

    lqn::LqnModel<T> model_;
    LqnsOptions opt_;
    lqn::LqnStruct<T> sn_;
    LqnsRawAvg raw_;
    std::map<std::pair<int, std::string>, std::size_t> byname_;
    std::map<std::string, std::size_t> bycallname_;
    double runtime_ = 0.0;
    bool solved_ = false;
};

}  // namespace lqns
}  // namespace line

#endif  // LINE_SOLVERS_WRAPPERS_LQNS_SOLVER_LQNS_H
