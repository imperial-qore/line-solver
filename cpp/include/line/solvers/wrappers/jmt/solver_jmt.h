/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_WRAPPERS_JMT_SOLVER_JMT_H
#define LINE_SOLVERS_WRAPPERS_JMT_SOLVER_JMT_H

/**
 * Port of `SolverJMT`, the Java Modelling Tools client.
 *
 * WHAT A JMT SOLVE IS. The model is written out as a JMT document -- a
 * `.jsimg` for the discrete-event engine, a `.jmva` for the analytical one --
 * `jmt.commandline.Jmt` is run on it, and the result document JMT leaves beside
 * the model is parsed back into the same `AvgResult` every other solver in this
 * port returns. The translation is `io/jmt_writer.h` and the SHARED
 * `io/jmva_writer.h`, which SolverQNS writes through as well; this file is the
 * dispatch, the parse and the metric mapping.
 *
 * THREE WAYS TO REACH JMT, in the order `jmtRun.m` tries them:
 *   1. `options.rest_url`, a JMT REST server (the imperialqore/jmt-rest image).
 *      Nothing runs locally.
 *   2. a local JVM plus `common/JMT.jar`, the default.
 *   3. no JVM, but a usable Docker daemon: the same image, run as a container.
 * Every one of them leaves the result at `<model>-result.jsim` or
 * `<model>-result.jmva`, so the parsers do not know which one ran. That
 * contract is what makes the three interchangeable, and it is why the Docker
 * arm copies the container's output back beside the model.
 *
 * CONSENT IS NOT ASSUMED FOR THE PULL. The Docker arm is reached only when
 * `LINE_JMT_DOCKER` opts in: this port has no terminal to ask at -- the
 * reference prompts, and refuses in a `-batch` session -- so an absent variable
 * is a refusal, exactly as an empty answer is there. Downloading 50 MB because
 * a solver was called is not a decision a library may take.
 *
 * WHAT IS NOT PORTED, and why it is refused rather than approximated:
 *   getProbAggr / getProbSysAggr   Both weigh the simulated trajectory against
 *                                  `sn.state{isf}`, the model's CURRENT state,
 *                                  which `qn::NetworkStruct` does not carry --
 *                                  the same gap that stops `SolverLDES`'s
 *                                  `getProb` in this port.
 *   the `replication` method       It averages `sampleSysAggr` over `iter_max`
 *                                  seeds, and `sampleSysAggr` reads the JMT
 *                                  arrival/departure LOG files, which requires
 *                                  a Logger on every station of the model; the
 *                                  transient arm is available through
 *                                  `jmt_sample_sys_aggr` and is composed there.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

#include "line/io/docker_image.h"
#include "line/io/jmva_writer.h"
#include "line/io/jmt_writer.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"
#include "line/util/http.h"
#include "line/util/subprocess.h"
#include "line/util/tempdir.h"
#include "line/util/xml.h"

namespace line {
namespace jmt {

/** The default JMT REST/Docker image, MATLAB's `jmtDockerImage` candidate. */
static const char* JMT_DOCKER_IMAGE = "imperialqore/jmt-rest:latest";

/** The options of one JMT solve, `SolverOptions('JMT')` restricted to what is read. */
struct JmtOptions {
    std::string method = "default";  ///< default | jsim | jmva | jmva.<alg>
    double samples = 10000.0;        ///< samples per measure; raised to 5000 below
    long seed = 23000;
    bool keep = false;               ///< keep the scratch directory after the solve
    double confint = 0.99;           ///< 0 disables the confidence interval columns
    double max_simulated_time = std::numeric_limits<double>::infinity();
    std::string rest_url;            ///< a JMT REST server; empty selects the local JVM
    std::string container;           ///< `options.config.container`, an image override
    int timeout = 3600;              ///< seconds, for the REST and subprocess arms
    bool verbose = false;
    /** `options.iter_max`: the replication count of `method = "replication"`. */
    int iter_max = 10;
};

namespace detail {

inline bool is_file(const std::string& p) {
    struct stat st;
    return !p.empty() && ::stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

inline bool is_exec(const std::string& p) {
    struct stat st;
    return !p.empty() && ::stat(p.c_str(), &st) == 0 && ::access(p.c_str(), X_OK) == 0;
}

inline std::string env_or_empty(const char* name) {
    const char* v = std::getenv(name);
    return v != nullptr ? std::string(v) : std::string();
}

inline std::string exe_dir() {
    char buf[4096];
    const ssize_t n = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n <= 0) return std::string();
    std::string p(buf, static_cast<std::size_t>(n));
    const std::size_t s = p.find_last_of('/');
    return s == std::string::npos ? std::string() : p.substr(0, s);
}

inline std::string cwd() {
    char buf[4096];
    return ::getcwd(buf, sizeof(buf)) != nullptr ? std::string(buf) : std::string();
}

/**
 * Port of `jmtGetPath`, minus the download.
 *
 * The reference fetches JMT.jar from SourceForge when it is missing. THIS PORT
 * DOES NOT: a solver call is not consent to a 50 MB download, and the reference
 * only gets away with it because it prints and prompts. `jmt_jar_path` returns
 * empty instead, and the dispatch then says what to install and where.
 *
 * `$LINE_JMT_DIR` overrides the search, for a jar kept outside the tree.
 */
inline std::string jmt_jar_path() {
    const std::string forced = env_or_empty("LINE_JMT_DIR");
    if (!forced.empty()) return is_file(forced + "/JMT.jar") ? forced : std::string();
    const std::string roots[2] = {exe_dir(), cwd()};
    for (int r = 0; r < 2; ++r) {
        if (roots[r].empty()) continue;
        std::string dir = roots[r];
        for (int k = 0; k < 7; ++k) {
            if (is_file(dir + "/JMT.jar")) return dir;
            if (is_file(dir + "/common/JMT.jar")) return dir + "/common";
            dir += "/..";
        }
    }
    return std::string();
}

/** `$LINE_JAVA`, then `$JAVA_HOME/bin/java`, then PATH; empty when absent. */
inline std::string find_java() {
    const std::string forced = env_or_empty("LINE_JAVA");
    if (!forced.empty()) return is_exec(forced) ? forced : std::string();
    const std::string home = env_or_empty("JAVA_HOME");
    if (!home.empty() && is_exec(home + "/bin/java")) return home + "/bin/java";
    const std::string path = env_or_empty("PATH");
    std::size_t b = 0;
    while (b <= path.size()) {
        const std::size_t e = path.find(':', b);
        const std::string dir =
            path.substr(b, e == std::string::npos ? std::string::npos : e - b);
        if (!dir.empty() && is_exec(dir + "/java")) return dir + "/java";
        if (e == std::string::npos) break;
        b = e + 1;
    }
    return std::string();
}

/** `LINE_JMT_DOCKER` in the affirmative; anything else, including unset, is no. */
inline bool docker_consented() {
    const std::string v = env_or_empty("LINE_JMT_DOCKER");
    return v == "1" || v == "true" || v == "TRUE" || v == "yes" || v == "y" || v == "Y";
}

inline std::string read_file(const std::string& path) {
    std::ifstream f(path.c_str(), std::ios::binary);
    if (!f) throw InputError("SolverJMT: cannot read '" + path + "'");
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

inline void write_file(const std::string& path, const std::string& text) {
    std::ofstream f(path.c_str(), std::ios::binary);
    if (!f) throw InputError("SolverJMT: cannot write '" + path + "'");
    f.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!f) throw InputError("SolverJMT: write failed on '" + path + "'");
}

/** JSON string escaping, for the REST request body. */
inline std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 16);
    for (std::size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                } else {
                    out.push_back(s[i]);
                }
        }
    }
    return out;
}

/**
 * The value of a top-level JSON string field, unescaped.
 *
 * The REST response is a small, flat envelope -- `status`, `error`,
 * `raw_output.result_xml`, `raw_output.stdout` -- and pulling four strings out
 * of it does not warrant making this header depend on the JSON library that
 * only the LDES client uses. A malformed response yields an empty string, and
 * the caller then reports that the server sent no result document.
 */
inline std::string json_string_field(const std::string& body, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    std::size_t p = body.find(needle);
    if (p == std::string::npos) return std::string();
    p = body.find(':', p + needle.size());
    if (p == std::string::npos) return std::string();
    ++p;
    while (p < body.size() && (body[p] == ' ' || body[p] == '\t' || body[p] == '\n')) ++p;
    if (p >= body.size() || body[p] != '"') return std::string();
    ++p;
    std::string out;
    while (p < body.size() && body[p] != '"') {
        if (body[p] == '\\' && p + 1 < body.size()) {
            ++p;
            switch (body[p]) {
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                case 'u': {
                    if (p + 4 < body.size()) {
                        const long cp = std::strtol(body.substr(p + 1, 4).c_str(), nullptr, 16);
                        if (cp < 0x80) out.push_back(static_cast<char>(cp));
                        p += 4;
                    }
                    break;
                }
                default: out.push_back(body[p]);
            }
        } else {
            out.push_back(body[p]);
        }
        ++p;
    }
    return out;
}

}  // namespace detail

template <class T>
struct JmtResult;

inline util::ProcResult jmt_run_docker(const std::string& image, const std::string& mode,
                                       const std::string& model_path, long seed,
                                       const JmtOptions& opt);

/** The suffix `jmt.commandline.Jmt` appends to the model path, per mode. */
inline const char* jmt_result_ext(const std::string& mode) {
    if (mode == "sim") return "jsim";
    if (mode == "mva") return "jmva";
    throw InputError("SolverJMT: unknown JMT analysis mode '" + mode + "'");
}

/** True when a local JVM and `common/JMT.jar` are both present. */
inline bool jmt_available() {
    return !detail::jmt_jar_path().empty() && !detail::find_java().empty();
}

/**
 * Port of `jmtSolveRest`: POST the model document, write the result document
 * back beside the model so the local parsers are unaffected.
 */
inline util::ProcResult jmt_solve_rest(const std::string& rest_url, const std::string& mode,
                                       const std::string& model_path, long seed,
                                       const JmtOptions& opt) {
    std::string url = rest_url;
    while (!url.empty() && url[url.size() - 1] == '/') url.erase(url.size() - 1);
    // A caller may hand over the full route; anything else gets the default one.
    if (url.size() < 8 || url.substr(url.size() - 8) != "/solve/" + mode.substr(0, 1))
        if (url.find("/api/v") == std::string::npos) url += "/api/v1/solve/" + mode;

    std::string body = "{\"model\":{\"content\":\"" +
                       detail::json_escape(detail::read_file(model_path)) +
                       "\",\"base64\":false}";
    // JMVA takes its algorithm and tolerance from the model document, so the
    // seed is meaningful only on the simulation route; the server's option
    // allow-list rejects it on /solve/mva.
    if (mode == "sim") body += ",\"options\":{\"seed\":" + std::to_string(seed) + "}";
    body += "}";

    const http::Response resp = http::post_json(url, body, opt.timeout * 1000);
    if (resp.status != 200)
        throw InputError("SolverJMT: the JMT REST server at " + url + " answered HTTP " +
                         std::to_string(resp.status));
    const std::string status = detail::json_string_field(resp.body, "status");
    if (status != "completed") {
        const std::string err = detail::json_string_field(resp.body, "error");
        throw InputError("SolverJMT: the JMT REST solve failed: " +
                         (err.empty() ? std::string("unspecified error") : err));
    }
    const std::string xml = detail::json_string_field(resp.body, "result_xml");
    if (xml.empty())
        throw InputError(
            "SolverJMT: the JMT REST response carries no result document. The server was asked "
            "to include the raw output; check that include_raw_output is not disabled");
    detail::write_file(model_path + "-result." + jmt_result_ext(mode), xml);
    util::ProcResult r;
    r.exitCode = 0;
    r.out = detail::json_string_field(resp.body, "stdout");
    return r;
}

/**
 * Port of `jmtRun`: one batch analysis, leaving the result where the JMT CLI
 * itself would leave it.
 *
 * @param mode "sim" or "mva"
 * @param model_path the document JMT is to read; the result lands beside it
 */
inline util::ProcResult jmt_run(const std::string& mode, const std::string& model_path, long seed,
                                const JmtOptions& opt) {
    if (!opt.rest_url.empty()) return jmt_solve_rest(opt.rest_url, mode, model_path, seed, opt);

    const std::string jar_dir = detail::jmt_jar_path();
    const std::string java = detail::find_java();
    if (!jar_dir.empty() && !java.empty()) {
        std::vector<std::string> argv;
        argv.push_back(java);
        argv.push_back("-cp");
        argv.push_back(jar_dir + "/JMT.jar");
        argv.push_back("jmt.commandline.Jmt");
        argv.push_back(mode);
        argv.push_back(model_path);
        argv.push_back("-seed");
        argv.push_back(std::to_string(seed));
        argv.push_back("--illegal-access=permit");
        const util::ProcResult r = util::capture(argv, opt.timeout, true);
        // A KILLED RUN IS NOT A RESULT. Without this the caller falls through to
        // "no result document", which reads as a solver that refused the model
        // rather than one whose wall-clock budget expired -- and a jsim that is
        // merely SLOW on a loaded host is exactly the case that hits it. Same
        // rule the LQNS, QNS and LDES wrappers already apply to util::capture.
        if (r.timedOut)
            throw NumericError("SolverJMT: the JMT run did not finish within " +
                               std::to_string(opt.timeout) + "s and was killed");
        return r;
    }

    // The Docker arm. The image is used only when it is ALREADY PRESENT or the
    // caller has opted into the pull; see the header note on consent.
    std::string image = opt.container.empty() ? detail::env_or_empty("LINE_JMT_IMAGE")
                                              : opt.container;
    if (image.empty()) image = JMT_DOCKER_IMAGE;
    if (io::docker_daemon_available()) {
        bool have = io::docker_has_local_image(image);
        if (!have && detail::docker_consented() && io::docker_has_storage_for(image))
            have = io::docker_pull(image);
        if (have) return jmt_run_docker(image, mode, model_path, seed, opt);
    }

    throw UnsupportedError(
        "SolverJMT: no way to reach JMT. Install a Java runtime and put JMT.jar in common/ "
        "(https://line-solver.sourceforge.net/latest/JMT.jar), or point options.rest_url at a "
        "JMT REST server, or set LINE_JMT_DOCKER=1 to allow the " +
        image + " image to be pulled and used");
}

/**
 * Run the analysis inside the JMT container and copy the result back beside the
 * model, so the caller sees the layout a local JVM would have produced.
 *
 * The model is STAGED under a fresh directory and the container is given that
 * directory: in `mva` mode JMT rewrites the model file itself, and a confined
 * Docker cannot bind-mount the system temp directory the model normally lives
 * in. The container runs as the calling user so the files it writes are not
 * left owned by root.
 */
inline util::ProcResult jmt_run_docker(const std::string& image, const std::string& mode,
                                       const std::string& model_path, long seed,
                                       const JmtOptions& opt) {
    util::TempDir work("jmt-docker", true);
    const std::size_t slash = model_path.find_last_of('/');
    const std::string base =
        slash == std::string::npos ? model_path : model_path.substr(slash + 1);
    const std::string staged = work.path() + "/" + base;
    detail::write_file(staged, detail::read_file(model_path));

    const std::string uid = std::to_string(static_cast<long>(::getuid()));
    const std::string gid = std::to_string(static_cast<long>(::getgid()));
    std::vector<std::string> argv;
    argv.push_back("docker");
    argv.push_back("run");
    argv.push_back("--rm");
    argv.push_back("--user");
    argv.push_back(uid + ":" + gid);
    argv.push_back("-v");
    argv.push_back(work.path() + ":" + work.path());
    argv.push_back("-w");
    argv.push_back(work.path());
    argv.push_back(image);
    argv.push_back(mode);
    argv.push_back(base);
    argv.push_back("-seed");
    argv.push_back(std::to_string(seed));
    const util::ProcResult r = util::capture(argv, opt.timeout, true);
    if (r.timedOut)
        throw NumericError("SolverJMT: the JMT container run did not finish within " +
                           std::to_string(opt.timeout) + "s and was killed");

    const std::string produced = staged + "-result." + jmt_result_ext(mode);
    if (detail::is_file(produced))
        detail::write_file(model_path + "-result." + jmt_result_ext(mode),
                           detail::read_file(produced));
    return r;
}

/** One `<measure>` of a JMT result document, by its attributes. */
struct JmtMeasure {
    std::string measure_type, station, job_class, node_type;
    double mean = 0.0, lower = 0.0, upper = 0.0;
    double analyzed_samples = 0.0;
    bool successful = false;
    bool has_bounds = false;
};

/**
 * Port of `getResultsJSIM`: every `<measure>` of the result document.
 *
 * A MISSING FILE IS A FAILED SIMULATION, not an empty result: JMT writes the
 * document even when every measure is unsuccessful, so its absence means the
 * run itself did not complete, and returning zeros would report a solved model.
 */
inline std::vector<JmtMeasure> jmt_parse_measures(const std::string& result_path,
                                                  const std::string& command_output) {
    if (!detail::is_file(result_path)) {
        std::string msg =
            "SolverJMT: JMT did not output a result file, the simulation has likely failed.";
        if (!command_output.empty()) msg += " JMT output: " + command_output;
        throw NumericError(msg);
    }
    const std::unique_ptr<xml::Element> doc = xml::parse_file(result_path);
    std::vector<JmtMeasure> out;
    const std::vector<const xml::Element*> ms = doc->by_tag("measure");
    for (std::size_t i = 0; i < ms.size(); ++i) {
        const xml::Element* e = ms[i];
        JmtMeasure m;
        m.measure_type = e->attr("measureType");
        m.station = e->attr("station");
        m.job_class = e->attr("class");
        m.node_type = e->attr("nodeType");
        m.mean = std::strtod(e->attr("meanValue").c_str(), nullptr);
        m.analyzed_samples = std::strtod(e->attr("analyzedSamples").c_str(), nullptr);
        m.successful = e->attr("successful") == "true";
        if (e->has_attr("lowerLimit") && e->has_attr("upperLimit")) {
            m.lower = std::strtod(e->attr("lowerLimit").c_str(), nullptr);
            m.upper = std::strtod(e->attr("upperLimit").c_str(), nullptr);
            m.has_bounds = true;
        }
        out.push_back(m);
    }
    return out;
}

/**
 * The result of a JMT solve: the shared `AvgResult` plus what only JMT reports.
 *
 * The FCR rows extend the metric matrices past `nstations`, which is the
 * reference's layout (`getResults.m` allocates `nstations + nregions` rows) and
 * is what makes a region's measures reachable through the same table as a
 * station's.
 */
template <class T>
struct JmtResult {
    mva::AvgResult<T> avg;
    Matrix<T> QCI, UCI, RCI, TCI, ACI;   ///< confidence half-widths, empty when disabled
    Matrix<T> Weight, MemOcc;            ///< FCR-only: weighted and memory occupation
    Matrix<T> TNfcr, DropRateNfcr;       ///< (nregions x nclasses) carried and lost rate
    /** Per Cache node (1-based node index), the per-class hit probability. */
    std::map<std::size_t, std::vector<T>> cache_hit_prob;
    /**
     * `result.Prob.logNormConstAggr`, the log normalizing constant JMVA reports.
     *
     * NaN on the simulation path and on any JMVA algorithm that does not
     * compute one -- an AMVA approximation has no G to report -- rather than
     * zero, which is a legitimate value of a log constant.
     */
    T log_norm_const = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
};

/**
 * Port of `getResults.m`: the measures mapped onto the metric matrices.
 *
 * THE RECURRENCE TEST IS WHY `analyzedSamples` IS READ. For a CLOSED class JMT
 * may terminate a measure before it has seen more samples than there are jobs
 * in the chain, and such a measure is not an estimate of anything; the
 * reference zeroes it rather than reporting it, and so does this. An OPEN class
 * has no such bound and is taken as reported.
 */
template <class T>
JmtResult<T> jmt_map_measures(const qn::NetworkStruct<T>& sn,
                              const std::vector<JmtMeasure>& measures, bool confint) {
    const std::size_t M = sn.nstations, K = sn.nclasses, F = sn.regions.size();
    const T zero = num_traits<T>::from_int(0);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    JmtResult<T> res;
    res.avg.QN = Matrix<T>(M + F, K, zero);
    res.avg.UN = Matrix<T>(M + F, K, zero);
    res.avg.RN = Matrix<T>(M + F, K, zero);
    res.avg.TN = Matrix<T>(M + F, K, zero);
    res.avg.AN = Matrix<T>(M + F, K, zero);
    res.avg.WN = Matrix<T>(M + F, K, zero);
    res.Weight = Matrix<T>(M + F, K, num_traits<T>::from_double(nan));
    res.MemOcc = Matrix<T>(M + F, K, num_traits<T>::from_double(nan));
    if (confint) {
        res.QCI = Matrix<T>(M, K, zero);
        res.UCI = Matrix<T>(M, K, zero);
        res.RCI = Matrix<T>(M, K, zero);
        res.TCI = Matrix<T>(M, K, zero);
        res.ACI = Matrix<T>(M, K, zero);
    }
    // JMT reports no utilization and no arrival rate for a region.
    for (std::size_t f = 0; f < F; ++f)
        for (std::size_t r = 0; r < K; ++r) {
            res.avg.UN(M + f, r) = num_traits<T>::from_double(nan);
            res.avg.AN(M + f, r) = num_traits<T>::from_double(nan);
        }

    std::map<std::string, std::size_t> station_of, class_of;
    for (std::size_t i = 1; i <= sn.nodes.size(); ++i) station_of[sn.nodes[i - 1].name] = i;
    for (std::size_t r = 1; r <= K; ++r) class_of[sn.classes[r - 1].name] = r;

    // The chain population each closed class's recurrence test compares against.
    std::vector<double> chainpop(K, 0.0);
    for (std::size_t c = 0; c < sn.inchain.size(); ++c) {
        double tot = 0.0;
        for (std::size_t r : sn.inchain[c]) {
            const double n = sn.classes[r - 1].population;
            if (std::isfinite(n)) tot += n;
        }
        for (std::size_t r : sn.inchain[c]) chainpop[r - 1] = tot;
    }

    for (std::size_t i = 0; i < measures.size(); ++i) {
        const JmtMeasure& m = measures[i];
        const double half = m.has_bounds ? (m.upper - m.lower) / 2.0 : 0.0;

        // A region measure is identified by its nodeType, or by the name the
        // writer gave it -- JMT does not always echo the nodeType back.
        const bool is_fcr =
            m.node_type == "region" || m.station.compare(0, 8, "FCRegion") == 0;
        if (is_fcr && m.station.compare(0, 8, "FCRegion") == 0) {
            const long f = std::strtol(m.station.substr(8).c_str(), nullptr, 10);
            if (f < 1 || static_cast<std::size_t>(f) > F) continue;
            const std::size_t row = M + static_cast<std::size_t>(f) - 1;
            // A region measure is AGGREGATE, not per class. The reference
            // spreads the extensive ones (customers, throughput, arrival rate,
            // weight, memory) evenly over the classes and repeats the intensive
            // ones (response and residence time) unchanged, so that summing a
            // column recovers the region total in both cases.
            const double per_class = m.mean / static_cast<double>(K);
            for (std::size_t r = 0; r < K; ++r) {
                if (m.measure_type == "Number of Customers")
                    res.avg.QN(row, r) = num_traits<T>::from_double(per_class);
                else if (m.measure_type == "Response Time")
                    res.avg.RN(row, r) = num_traits<T>::from_double(m.mean);
                else if (m.measure_type == "Residence Time")
                    res.avg.WN(row, r) = num_traits<T>::from_double(m.mean);
                else if (m.measure_type == "Throughput")
                    res.avg.TN(row, r) = num_traits<T>::from_double(per_class);
                else if (m.measure_type == "Arrival Rate")
                    res.avg.AN(row, r) = num_traits<T>::from_double(per_class);
                else if (m.measure_type == "FCR Capacity")
                    res.Weight(row, r) = num_traits<T>::from_double(per_class);
                else if (m.measure_type == "FCR Memory")
                    res.MemOcc(row, r) = num_traits<T>::from_double(per_class);
            }
            continue;
        }

        if (m.measure_type == "Cache Hit Rate") {
            const auto ni = station_of.find(m.station);
            const auto ci = class_of.find(m.job_class);
            if (ni == station_of.end() || ci == class_of.end()) continue;
            if (sn.nodes[ni->second - 1].nodetype != lang::NodeType::Cache) continue;
            const auto cp = sn.nodeparam.find(ni->second);
            if (cp == sn.nodeparam.end()) continue;
            std::vector<T>& hit = res.cache_hit_prob[ni->second];
            if (hit.empty()) hit.assign(K, zero);
            // The measure names the HIT class; the probability belongs to the
            // read class that switches into it.
            for (std::size_t r = 0; r < cp->second.hitclass.size(); ++r)
                if (cp->second.hitclass[r] == ci->second)
                    hit[r] = num_traits<T>::from_double(m.mean);
            continue;
        }

        const auto ni = station_of.find(m.station);
        const auto ci = class_of.find(m.job_class);
        if (ni == station_of.end() || ci == class_of.end()) continue;
        const std::size_t ist = sn.nodes[ni->second - 1].station;
        if (ist == 0) continue;
        const std::size_t r = ci->second;
        const bool open = !std::isfinite(sn.classes[r - 1].population);
        const bool recurrent = open || m.analyzed_samples > chainpop[r - 1];
        const T v = num_traits<T>::from_double(recurrent ? m.mean : 0.0);
        const T ci_v = num_traits<T>::from_double(recurrent ? half : 0.0);

        if (m.measure_type == "Number of Customers") {
            res.avg.QN(ist - 1, r - 1) = v;
            if (confint) res.QCI(ist - 1, r - 1) = ci_v;
        } else if (m.measure_type == "Utilization") {
            res.avg.UN(ist - 1, r - 1) = v;
            if (confint) res.UCI(ist - 1, r - 1) = ci_v;
        } else if (m.measure_type == "Response Time") {
            res.avg.RN(ist - 1, r - 1) = v;
            if (confint) res.RCI(ist - 1, r - 1) = ci_v;
        } else if (m.measure_type == "Throughput") {
            res.avg.TN(ist - 1, r - 1) = v;
            if (confint) res.TCI(ist - 1, r - 1) = ci_v;
        } else if (m.measure_type == "Arrival Rate") {
            res.avg.AN(ist - 1, r - 1) = v;
            if (confint) res.ACI(ist - 1, r - 1) = ci_v;
        }
        // `Residence Time` is deliberately neither requested nor read: JMT's
        // definition disagrees with LINE's on class-switching models, and the
        // residence time is recomputed below from the response time.
    }
    return res;
}

/**
 * The region loss table, port of the `sn.nregions > 0` tail of `getResults.m`.
 *
 * JMT EXPOSES NO REGION DROP MEASURE, and its region throughput is the CARRIED
 * (admitted) rate, so the offered rate is reconstructed by flow balance: the
 * rate at which stations outside the region route jobs across its boundary,
 * which the drop does not affect. The difference, clamped at zero, is the loss;
 * a region whose rule is not DROP loses nothing by construction and is zeroed.
 */
template <class T>
void jmt_region_losses(const qn::NetworkStruct<T>& sn, JmtResult<T>& res) {
    const std::size_t M = sn.nstations, K = sn.nclasses, F = sn.regions.size();
    if (F == 0) return;
    const T zero = num_traits<T>::from_int(0);
    res.TNfcr = Matrix<T>(F, K, zero);
    res.DropRateNfcr = Matrix<T>(F, K, zero);
    for (std::size_t f = 0; f < F; ++f)
        for (std::size_t r = 0; r < K; ++r) res.TNfcr(f, r) = res.avg.TN(M + f, r);

    for (std::size_t f = 0; f < F; ++f) {
        const typename qn::NetworkStruct<T>::Region& rg = sn.regions[f];
        for (std::size_t r = 1; r <= K; ++r) {
            double acc = 0.0;
            for (std::size_t ii = 1; ii <= M; ++ii) {
                if (!(ii <= rg.members.size() && rg.members[ii - 1])) continue;
                for (std::size_t j = 1; j <= M; ++j) {
                    if (j <= rg.members.size() && rg.members[j - 1]) continue;
                    for (std::size_t rp = 1; rp <= K; ++rp) {
                        const std::size_t a = (j - 1) * K + (rp - 1);
                        const std::size_t b = (ii - 1) * K + (r - 1);
                        if (sn.rt.rows() <= a || sn.rt.cols() <= b) continue;
                        const double w = num_traits<T>::to_double(sn.rt(a, b));
                        if (w != 0.0)
                            acc += num_traits<T>::to_double(res.avg.TN(j - 1, rp - 1)) * w;
                    }
                }
            }
            const double carried = num_traits<T>::to_double(res.TNfcr(f, r - 1));
            const bool drops = r <= rg.rule.size() && rg.rule[r - 1] == lang::DropStrategy::DROP;
            res.DropRateNfcr(f, r - 1) =
                num_traits<T>::from_double(drops ? std::max(0.0, acc - carried) : 0.0);
        }
    }
}

/**
 * Port of `getResultsJMVA.m`: the per-CHAIN answer spread back over the classes.
 *
 * JMVA answers per chain, and each measure is converted to a per-class one by
 * the class's share of its chain at that station: `alpha(i,k)` weighted by the
 * ratio of the class service time to the chain's. The conversions are the
 * reference's, measure by measure, and the two that are not a bare share are
 * the reason this cannot be a generic rescale:
 *   Utilization   is divided by the chain's visits at the REFERENCE station and,
 *                 at a multiserver station, scaled by min(N, c)/c.
 *   Residence time is JMVA's per-CHAIN residence and is converted to LINE's
 *                 response time per visit by dividing by the class visits.
 *
 * THE STATION IS RESOLVED BY NAME, not by the position of the `stationresults`
 * block. The reference indexes `sn.nservers`, `ST` and `sn.visits` with the
 * BLOCK index, but the blocks exclude the Source while those tables do not, so
 * on any open model every station's demand is read one row early. Reading the
 * `station` attribute is the same lookup on a closed model and the correct one
 * on an open model.
 */
template <class T>
JmtResult<T> jmt_parse_jmva(const qn::NetworkStruct<T>& sn, const std::string& result_path,
                            const std::string& command_output) {
    if (!detail::is_file(result_path)) {
        std::string msg =
            "SolverJMT: JMT did not output a result file, the analysis has likely failed.";
        if (!command_output.empty()) msg += " JMT output: " + command_output;
        throw NumericError(msg);
    }
    const std::size_t M = sn.nstations, K = sn.nclasses, C = sn.nchains;
    const T zero = num_traits<T>::from_int(0);
    JmtResult<T> res;
    res.avg.QN = Matrix<T>(M, K, zero);
    res.avg.UN = Matrix<T>(M, K, zero);
    res.avg.RN = Matrix<T>(M, K, zero);
    res.avg.TN = Matrix<T>(M, K, zero);
    res.avg.AN = Matrix<T>(M, K, zero);
    res.avg.WN = Matrix<T>(M, K, zero);

    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(sn);

    std::map<std::string, std::size_t> station_of;
    for (std::size_t i = 1; i <= M; ++i)
        station_of[sn.nodes[sn.station_to_node[i - 1] - 1].name] = i;

    const std::unique_ptr<xml::Element> doc = xml::parse_file(result_path);
    // The normalizing constant, when the algorithm reports one.
    const std::vector<const xml::Element*> nc = doc->by_tag("normconst");
    res.log_norm_const = nc.empty()
                             ? num_traits<T>::from_double(
                                   std::numeric_limits<double>::quiet_NaN())
                             : num_traits<T>::from_double(
                                   std::strtod(nc[0]->attr("logValue").c_str(), nullptr));

    // A Source reports its own arrival rate as its throughput and holds no
    // jobs; JMVA does not model it, so the row is filled from the struct.
    if (sn.sourceIdx != 0)
        for (std::size_t r = 0; r < K; ++r)
            if (!sn.disabled[sn.sourceIdx - 1][r])
                res.avg.TN(sn.sourceIdx - 1, r) = sn.rates(sn.sourceIdx - 1, r);

    const std::vector<const xml::Element*> blocks = doc->by_tag("stationresults");
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        const auto si = station_of.find(blocks[b]->attr("station"));
        if (si == station_of.end()) continue;
        const std::size_t i = si->second;
        // The divisor is the capacity `write_jmva` exported, max(nservers, max
        // lldscaling): a load-dependent station carries its c in the scaling and
        // leaves nservers at 1, so reading nservers alone reports U = c*E[busy]/c.
        double ns = sn.stations[i - 1].nservers;
        for (const T& s : sn.stations[i - 1].lldscaling)
            ns = std::max(ns, num_traits<T>::to_double(s));
        const std::vector<const xml::Element*> crs = blocks[b]->by_tag("classresults");
        for (std::size_t c = 1; c <= crs.size() && c <= C; ++c) {
            const std::vector<const xml::Element*> ms = crs[c - 1]->by_tag("measure");
            // A multiserver Queue is written as <ldstation servers="1">, and one
            // ldstation switches JMVA to its load-dependent algorithm, whose
            // Utilization is 1-p_i(0) at EVERY station, delay ones included. That
            // is a different random variable from LINE's E[busy servers], so no
            // rescaling recovers it; derive U from the chain throughput instead,
            // which both JMVA algorithms report alike.
            double chain_tput = std::numeric_limits<double>::quiet_NaN();
            for (std::size_t m = 0; m < ms.size(); ++m) {
                if (ms[m]->attr("measureType") == "Throughput") {
                    const std::string traw = ms[m]->attr("meanValue");
                    chain_tput = traw == "NaN" ? std::numeric_limits<double>::quiet_NaN()
                                               : std::strtod(traw.c_str(), nullptr);
                    break;
                }
            }
            for (std::size_t m = 0; m < ms.size(); ++m) {
                const std::string kind = ms[m]->attr("measureType");
                const std::string raw = ms[m]->attr("meanValue");
                const double val = raw == "NaN" ? std::numeric_limits<double>::quiet_NaN()
                                                : std::strtod(raw.c_str(), nullptr);
                const std::size_t refst = dem.refstatchain[c - 1];
                const double vchain_ref =
                    num_traits<T>::to_double(dem.Vchain(refst - 1, c - 1));
                const double stchain = num_traits<T>::to_double(dem.STchain(i - 1, c - 1));
                for (std::size_t k : sn.inchain[c - 1]) {
                    const double stk = num_traits<T>::to_double(dem.ST(i - 1, k - 1));
                    const double al = num_traits<T>::to_double(dem.alpha(i - 1, k - 1));
                    double v = val;
                    if (kind == "Utilization") {
                        if (vchain_ref == 0.0) continue;
                        v = stk * chain_tput / vchain_ref * al;
                        if (std::isfinite(ns)) v /= ns;
                        res.avg.UN(i - 1, k - 1) = num_traits<T>::from_double(v);
                    } else if (kind == "Throughput") {
                        res.avg.TN(i - 1, k - 1) = num_traits<T>::from_double(val * al);
                    } else if (kind == "Number of Customers") {
                        if (stchain == 0.0 || vchain_ref == 0.0) continue;
                        res.avg.QN(i - 1, k - 1) =
                            num_traits<T>::from_double(val * stk / stchain / vchain_ref * al);
                    } else if (kind == "Residence time" || kind == "Residence Time") {
                        if (stchain == 0.0 || vchain_ref == 0.0) continue;
                        double visits = 0.0;
                        if (c - 1 < sn.visits.size()) {
                            const std::size_t sfi = sn.stateful_of_station(i);
                            if (sfi != 0 && sn.visits[c - 1].rows() >= sfi)
                                visits = num_traits<T>::to_double(sn.visits[c - 1](sfi - 1, k - 1));
                        }
                        if (visits == 0.0) continue;
                        v = (val / visits) * stk / stchain / vchain_ref * al;
                        res.avg.RN(i - 1, k - 1) = num_traits<T>::from_double(v);
                    }
                }
            }
        }
    }
    return res;
}

/**
 * Port of `SolverJMT.listValidMethods`.
 *
 * `replication` is the reference's own TRANSIENT route: a single sample path is
 * not the transient mean E[N](t), there being no time-ergodicity at a fixed t,
 * so it averages `sampleSysAggr` over `iter_max` seeds. It is composed from
 * `jmt_sample_sys_aggr` in `jmt_logs.h` (`jmt_replication`) rather than from
 * `solver_jmt_run_analyzer`, which is why it is not an arm of the dispatch below.
 */
inline std::vector<std::string> jmt_list_valid_methods() {
    return {"default",   "jsim",       "replication", "jmva",     "jmva.amva", "jmva.mva",
            "jmva.recal", "jmva.comom", "jmva.chow", "jmva.bs",   "jmva.aql",
            "jmva.lin",   "jmva.dmlin"};
}

/**
 * The structural half of SolverJMT's method gate; empty when admissible.
 *
 * THE RULES A FLAT FEATURE SET CANNOT STATE, and each is named rather than
 * reported as a bare no:
 *   the transient horizon   'replication' averages `iter_max` sample paths over
 *                           [0,T]; a mean at an unstated horizon is not a
 *                           quantity, and the horizon is an OPTION rather than a
 *                           model feature.
 *   a single-server station the eight closed-form JMVA algorithms are
 *                           single-server only, which is why `write_jmva`
 *                           refuses the model rather than emitting an
 *                           <ldstation> the algorithm cannot read. A server
 *                           count has no feature name; the load-dependent half
 *                           of the same restriction DOES, and rides in
 *                           `jmt_feature_set` as an unset LoadDependence.
 *   the load-dependent shape JMT carries a scaling only as a server count, so
 *                           alpha(n) = min(n,c) is written exactly and nothing
 *                           else is.
 *
 * ONE PREDICATE, TWO CALLERS: `solver_jmt_run_analyzer` raises it, so a caller
 * naming the method by hand gets the sentence rather than a JMT stack trace, and
 * `autosolver::auto_family_refusal` returns it, so findSolver never offers the
 * pair. A second copy of any rule is how the gate and the run drift apart.
 */
template <class T>
std::string jmt_method_refusal(const qn::NetworkStruct<T>& sn, const std::string& method,
                               const JmtOptions& opt) {
    if (method == "replication" && !std::isfinite(opt.max_simulated_time))
        return "SolverJMT: the 'replication' method needs a finite timespan; a transient mean "
               "over an unstated horizon is not a quantity";
    if (qn::jmva_is_closed_only(method)) {
        double maxFinite = 0.0;
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            const double c = static_cast<double>(sn.stations[i].nservers);
            if (std::isfinite(c)) maxFinite = std::max(maxFinite, c);
        }
        if (maxFinite > 1.0) return "SolverJMT: " + method + " does not support multi-server stations";
    }
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const std::vector<T>& alpha = sn.stations[i].lldscaling;
        if (alpha.empty()) continue;
        double c = 0.0;
        for (const T& s : alpha) c = std::max(c, num_traits<T>::to_double(s));
        if (c == 1.0) continue;
        bool ok = (c >= 1.0) && (c == std::floor(c));
        for (std::size_t n = 0; ok && n < alpha.size(); ++n)
            ok = std::abs(num_traits<T>::to_double(alpha[n]) -
                          std::min(static_cast<double>(n + 1), c)) <= qn::GlobalConstants::Zero;
        if (!ok)
            return "SolverJMT: station '" + sn.stations[i].name +
                   "' uses a load-dependent scaling that is not the multiserver encoding alpha(n) "
                   "= min(n,c): JMT has no representation for it, since both the JSIM and the JMVA "
                   "writer carry the scaling as a server count, and the model would be solved at "
                   "the nominal service rate. Use SolverCTMC, SolverNC, SolverMVA or SolverSSA, "
                   "which read sn.lldscaling directly";
    }
    // A BINDING FINITE BUFFER, which neither engine can carry -- JSIM because no
    // JMT drop strategy reproduces LINE's blocking, JMVA because its document has
    // no capacity element at all. The verdict is the WRITER's own
    // (`io::JmtWriter::buffer_capacity_refusal`), asked here without letting it
    // raise, so the gate binds exactly where the writer binds.
    //
    // WHICH ENGINE IS ASKED ABOUT: a method that is neither engine's gets no
    // verdict rather than JSIM's rule applied to a run JSIM is not going to make.
    const bool is_jmva = method.compare(0, 4, "jmva") == 0;
    const bool is_jsim = method == "default" || method == "jsim" || method == "replication";
    if (is_jmva || is_jsim) {
        const std::string cap_reason = io::jmt_buffer_capacity_refusal(sn, is_jmva);
        if (!cap_reason.empty()) return cap_reason;
    }
    return std::string();
}

/**
 * Port of `@@SolverJMT/runAnalyzer.m`, the `jsim` and `jmva` arms.
 *
 * The sample count is RAISED to 5000 when the caller asks for less, as the
 * reference does: JMT terminates a measure on its own precision target and
 * needs that many samples per measure before it will report one, so a smaller
 * request silently yields unsuccessful measures rather than a faster run.
 */
template <class T>
JmtResult<T> solver_jmt_run_analyzer(const qn::NetworkStruct<T>& sn, const JmtOptions& opt_in) {
    JmtOptions opt = opt_in;
    if (opt.samples < 5000.0) opt.samples = 5000.0;
    const std::vector<std::string> valid = jmt_list_valid_methods();
    if (std::find(valid.begin(), valid.end(), opt.method) == valid.end())
        throw UnsupportedError("SolverJMT: unknown method '" + opt.method + "'");
    if (opt.method == "replication")
        throw UnsupportedError(
            "SolverJMT: 'replication' is the transient average over independent seeds and is "
            "composed from jmt_sample_sys_aggr, not from this steady-state dispatch; call "
            "jmt::jmt_replication (jmt_logs.h), which is what -a tran reaches");
    if (sn.has_immediate_feedback())
        throw UnsupportedError(
            "SolverJMT: JMT has no immediate-feedback semantics (a job re-entering service "
            "while HOLDING the server); use SolverCTMC, SolverSSA or SolverLDES");
    // The structural half of the gate, which `auto_family_refusal` also asks:
    // the single-server restriction of the closed-form JMVA algorithms and the
    // load-dependent shape JMT can carry. A second copy of either rule is how
    // the report and the run come to disagree.
    {
        const std::string refusal = jmt_method_refusal(sn, opt.method, opt);
        if (!refusal.empty()) throw UnsupportedError(refusal);
    }

    const bool is_mva = opt.method.compare(0, 4, "jmva") == 0;
    const std::string mode = is_mva ? "mva" : "sim";
    util::TempDir work(is_mva ? "jmva" : "jsim");
    if (opt.keep) work.keep();
    const std::string model_path = work.path() + "/model." + (is_mva ? "jmva" : "jsim");

    if (is_mva) {
        // `io/jmva_writer.h` is the shared port of `writeJMVA`: SolverQNS reads
        // the same grammar through `qnsolver`, and a second copy of the writer
        // here would be a second place for the chain aggregation to drift.
        io::write_jmva(sn, model_path, opt.method, static_cast<std::size_t>(opt.samples));
    } else {
        io::JmtWriteOptions wopt;
        wopt.file_name = "model";
        wopt.log_path = sn.log_path;
        wopt.seed = opt.seed;
        wopt.max_samples = opt.samples;
        wopt.max_simulated_time = opt.max_simulated_time;
        wopt.sim_conf_int = opt.confint > 0.0 ? opt.confint : 0.99;
        line::util::LineConsole::step("writing the JSIM model file");
        detail::write_file(model_path, io::jmt_write_jsim(sn, wopt));
        line::util::LineConsole::substep("model written to %s", model_path.c_str());
    }

    line::util::LineConsole::step("running the JMT simulation engine as a subprocess");
    const util::ProcResult run = jmt_run(mode, model_path, opt.seed, opt);
    const std::string result_path = model_path + "-result." + jmt_result_ext(mode);

    line::util::LineConsole::step("parsing the JMT result files");
    JmtResult<T> res;
    if (is_mva) {
        res = jmt_parse_jmva(sn, result_path, run.out);
    } else {
        res = jmt_map_measures(sn, jmt_parse_measures(result_path, run.out), opt.confint > 0.0);
        jmt_region_losses(sn, res);
    }
    // ResidT DERIVED on both arms as getAvg does (JSIM disagrees under class switching, JMVA per-chain); else column zero. Region rows keep JMT total.
    {
        const Matrix<T> WNst = mva::sn_get_residt_from_respt(sn, res.avg.RN);
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t r = 0; r < sn.nclasses; ++r) res.avg.WN(i, r) = WNst(i, r);
    }
    // A SOURCE HAS NO ARRIVALS TO ITSELF. JMT reports an "Arrival Rate" measure
    // at every station the model declares, the Source included, where the
    // reference reports 0: `getAvg.m:197-199` builds a zeroMask over the Source
    // stations and applies it to the ArvR column for EVERY solver. This port's
    // analytical arms already leave that entry at zero, so the JMT reader was
    // the one path that carried JMT's own measure through -- 0.79526 against 0
    // on fcr_mm1kdrop, 0.30214 and 0.4979 on fcr_constraints. Region rows are
    // past `nstations` and are left alone, as the residence time above leaves
    // them.
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].nodetype != qn::NodeType::Source) continue;
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            res.avg.AN(i, r) = num_traits<T>::from_int(0);
    }
    // A CELL WITH NO RESPONSE TIME HOLDS NO JOBS AND BUSIES NO SERVER.
    // `@@NetworkSolver/getAvg` builds `zeroMask = RN < 10*FineTol` and applies
    // it to the queue length and the utilization for EVERY solver; the
    // analytical arms of this port take it through `mva::filter_metric`, and
    // the JMT reader was the one path that did not. JMT keeps sampling a
    // starved class's utilization long after its response time has been
    // discarded as non-recurrent, so the cell arrives as a small positive
    // number where every other codebase reports 0: 0.0010617 against 0 on
    // prio_hol_closed (PSQueue, Class2), on an otherwise IDENTICAL sample path.
    // Region rows past `nstations` are left alone, as the two rules above are.
    {
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t r = 0; r < sn.nclasses; ++r)
                if (num_traits<T>::to_double(res.avg.RN(i, r)) <
                    10.0 * lang::GlobalConstants::FineTol) {
                    res.avg.QN(i, r) = num_traits<T>::from_int(0);
                    res.avg.UN(i, r) = num_traits<T>::from_int(0);
                }
    }
    res.avg.method = opt_in.method;
    res.avg.actualmethod = opt.method;
    res.avg.iter = 1;
    return res;
}

}  // namespace jmt
}  // namespace line

#endif  // LINE_SOLVERS_WRAPPERS_JMT_SOLVER_JMT_H
