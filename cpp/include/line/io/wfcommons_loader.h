/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_WFCOMMONS_LOADER_H
#define LINE_IO_WFCOMMONS_LOADER_H

/**
 * @file
 * @ingroup line_io
 * Port of `matlab/src/io/WfCommonsLoader.m` (and `jline.io.WfCommonsLoader` /
 * `WfCommonsOptions`): a WfCommons workflow trace
 * (https://github.com/wfcommons/workflow-schema, schema 1.3 to 1.5) read into a
 * `workflow::Workflow<T>`.
 *
 * WHAT IS BUILT. One activity per task, named by the task's `id` (its `name`
 * when there is no `id`, `task_<i>` when neither), with a host demand fitted to
 * the task's `runtimeInSeconds` from the execution data: embedded in the task
 * for schema 1.4 and earlier, under `workflow.execution.tasks` for 1.5+. The
 * edges are the union of the `children` (task -> child) and `parents` (parent ->
 * task) lists, each pair once and children-derived edges first, so a trace that
 * lists parents only (Montage dss) gives the same graph. They become precedences
 * in two passes, exactly as the reference:
 * a task with several children whose children ALL lead directly to one common
 * descendant becomes an AND-fork / AND-join pair, and every edge left over
 * becomes a serial precedence.
 *
 * CHILDREN THAT NAME A TASK BY `name` WHEN THE TASKS CARRY AN `id`. The Pegasus
 * traces of schema 1.4 do this: task `mProject_ID0000001` has
 * `"id": "ID0000001"` and lists its children as `"mBackground_ID0000013"`. The
 * reference keys the lookup on the id alone, so EVERY edge of such a trace is
 * dropped without a word and the workflow reduces to its tasks in series (the
 * JAR on montage-chameleon-2mass-005d-001.json: 58 activities, 0 precedences).
 * `WfCommonsOptions::resolve_children_by_name` (default true) looks a child that
 * matches no id up among the task names; a child matching neither is dropped,
 * as in the reference. Parent references resolve by the same rule, and one
 * warning counts the unresolved references of both kinds. Set it false for the
 * id-only lookup. On a document whose children are ids, which is every 1.5
 * document, the two agree.
 *
 * `loadFromUrl` fetches http:// over `line/util/http.h`, the port's own HTTP/1.1
 * client, which has no TLS and follows no redirects. An https:// URL, which is
 * how the WfCommons repositories are served and what MATLAB's `webread` takes,
 * is handed to `curl` (else `wget`) on PATH, as SolverJMT fetches JMT.jar: no
 * TLS library becomes a build dependency, and a host with neither tool is told
 * so by name. The fetch stays on https across redirects.
 *
 * WARNINGS go to stderr with the `[LINE] Warning:` prefix the io layer uses;
 * the port has no `line_warning` channel.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "json.hpp"
#include "line/lang/dist_fitters.h"
#include "line/lang/lang_types.h"
#include "line/lang/workflow/workflow.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/http.h"
#include "line/util/subprocess.h"
#include "line/util/tempdir.h"

namespace line {
namespace io {

/** `options.distributionType`, the JAR's `WfCommonsOptions.DistributionType`. */
enum class WfDistributionType { EXP, DET, APH, HYPEREXP };

/**
 * `lower(options.distributionType)` as the reference switches on it: `exp`, `det`,
 * `aph`, `hyperexp`. Anything else is EXP, the reference's `otherwise` branch.
 */
inline WfDistributionType wf_distribution_type_from_string(const std::string& s) {
    std::string l(s);
    for (char& c : l) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (l == "det") return WfDistributionType::DET;
    if (l == "aph") return WfDistributionType::APH;
    if (l == "hyperexp") return WfDistributionType::HYPEREXP;
    return WfDistributionType::EXP;
}

/** The loader options; the defaults are the reference's `parseOptions`. */
struct WfCommonsOptions {
    WfDistributionType distribution_type = WfDistributionType::EXP;  ///< `distributionType`
    double default_scv = 1.0;         ///< `defaultSCV`, for APH and HyperExp
    double default_runtime = 1.0;     ///< `defaultRuntime`, when a task has no runtime
    bool use_execution_data = true;   ///< `useExecutionData`
    bool store_metadata = true;       ///< `storeMetadata`, onto `WorkflowActivity::metadata()`
    /** Resolve a child matching no task id by task name (see the file comment). */
    bool resolve_children_by_name = true;
};

namespace wfcommons_detail {

using json = nlohmann::json;

/** `WfCommonsLoader.SUPPORTED_SCHEMA_VERSIONS`. */
inline bool supported_version(const std::string& v) {
    return v == "1.3" || v == "1.4" || v == "1.5";
}

inline bool has(const json& o, const char* key) { return o.is_object() && o.contains(key); }

/** `~isempty(x)` on a decoded JSON value. */
inline bool nonempty(const json& v) {
    if (v.is_null()) return false;
    if (v.is_array() || v.is_object() || v.is_string()) return !v.empty();
    return true;
}

/** A JSON scalar as the text `strcmp` compares: strings verbatim, numbers as written. */
inline std::string as_text(const json& v) { return v.is_string() ? v.get<std::string>() : v.dump(); }

/** Port of `validateSchema`. */
inline void validate_schema(const json& data) {
    if (!has(data, "schemaVersion"))
        throw InputError("WfCommonsLoader: Missing schemaVersion field in WfCommons JSON.");
    const std::string version = as_text(data.at("schemaVersion"));
    if (!supported_version(version))
        std::cerr << "[LINE] Warning: WfCommonsLoader: Schema version " << version
                  << " may not be fully supported." << std::endl;
    if (!has(data, "workflow"))
        throw InputError("WfCommonsLoader: Missing workflow field in WfCommons JSON.");
    const json& wf = data.at("workflow");
    bool has_tasks = false;
    if (has(wf, "specification")) {
        const json& spec = wf.at("specification");
        has_tasks = has(spec, "tasks") && nonempty(spec.at("tasks"));
    } else if (has(wf, "tasks") && nonempty(wf.at("tasks"))) {
        has_tasks = true;
    }
    if (!has_tasks) throw InputError("WfCommonsLoader: Workflow must have at least one task.");
}

/** `[~, name, ~] = fileparts(s)`. */
inline std::string file_stem(const std::string& s) {
    const std::size_t slash = s.find_last_of("/\\");
    std::string f = slash == std::string::npos ? s : s.substr(slash + 1);
    const std::size_t dot = f.find_last_of('.');
    if (dot != std::string::npos) f = f.substr(0, dot);
    return f;
}

/** Port of `extractName`: the document's name, else the default's stem, sanitized. */
inline std::string extract_name(const json& data, const std::string& default_name) {
    std::string name;
    if (has(data, "name") && data.at("name").is_string() && !data.at("name").get<std::string>().empty())
        name = data.at("name").get<std::string>();
    else
        name = file_stem(default_name);
    for (char& c : name)
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_')) c = '_';
    if (name.empty()) name = "Workflow";
    return name;
}

/** Port of `getTaskId`. */
inline std::string task_id(const json& task, std::size_t idx1) {
    if (has(task, "id")) return as_text(task.at("id"));
    if (has(task, "name")) return as_text(task.at("name"));
    return "task_" + std::to_string(idx1);
}

/** Port of `fitDistribution`. */
template <class T>
lang::Distrib<T> fit_distribution(double runtime, const WfCommonsOptions& opt) {
    typedef lang::Distrib<T> D;
    if (runtime <= lang::GlobalConstants::FineTol) return D::immediate();
    const T m = num_traits<T>::from_double(runtime);
    switch (opt.distribution_type) {
        case WfDistributionType::DET:
            return D::det(m);
        case WfDistributionType::APH:
            return lang::aph_fit_mean_scv(m, num_traits<T>::from_double(opt.default_scv));
        case WfDistributionType::HYPEREXP:
            if (opt.default_scv > 1.0)
                return lang::hyperexp_fit_mean_scv(m, num_traits<T>::from_double(opt.default_scv));
            return D::exp_mean(m);
        case WfDistributionType::EXP:
        default:
            return D::exp_mean(m);
    }
}

/** Port of `extractMetadata`; values are the JSON text of each field. */
inline std::map<std::string, std::string> extract_metadata(const json& task, const std::string& id,
                                                           const json* exec) {
    std::map<std::string, std::string> md;
    md["taskId"] = json(id).dump();
    static const char* kTask[] = {"name", "inputFiles", "outputFiles"};
    for (const char* f : kTask)
        if (has(task, f)) md[f] = task.at(f).dump();
    if (exec) {
        static const char* kExec[] = {"executedAt",  "command",       "coreCount",   "avgCPU",
                                      "readBytes",   "writtenBytes",  "memoryInBytes",
                                      "energyInKWh", "avgPowerInW",   "priority",    "machines"};
        for (const char* f : kExec)
            if (has(*exec, f)) md[f] = exec->at(f).dump();
    }
    return md;
}

/** Port of `getReachableNodes`: BFS order, the start excluded. */
inline std::vector<std::size_t> reachable_nodes(std::size_t start,
                                                const std::vector<std::vector<std::size_t>>& adj) {
    std::vector<bool> visited(adj.size(), false);
    std::deque<std::size_t> q(1, start);
    visited[start] = true;
    std::vector<std::size_t> out;
    while (!q.empty()) {
        const std::size_t cur = q.front();
        q.pop_front();
        for (std::size_t nx : adj[cur])
            if (!visited[nx]) {
                visited[nx] = true;
                q.push_back(nx);
                out.push_back(nx);
            }
    }
    return out;
}

/**
 * Port of `findCommonJoin`. `intersect` returns its result SORTED, so the first
 * qualifying join is the one with the smallest task index.
 *
 * @return the join's 0-based index, or `adj.size()` when there is none
 */
inline std::size_t find_common_join(const std::vector<std::size_t>& children,
                                    const std::vector<std::vector<std::size_t>>& adj,
                                    const std::vector<std::size_t>& indeg) {
    const std::size_t none = adj.size();
    if (children.empty()) return none;
    std::vector<std::size_t> common = reachable_nodes(children[0], adj);
    std::sort(common.begin(), common.end());
    common.erase(std::unique(common.begin(), common.end()), common.end());
    for (std::size_t c = 1; c < children.size(); ++c) {
        std::vector<std::size_t> r = reachable_nodes(children[c], adj);
        std::sort(r.begin(), r.end());
        std::vector<std::size_t> both;
        std::set_intersection(common.begin(), common.end(), r.begin(), r.end(),
                              std::back_inserter(both));
        common.swap(both);
    }
    for (std::size_t node : common) {
        if (indeg[node] < children.size()) continue;
        bool all_direct = true;
        for (std::size_t ch : children)
            if (std::find(adj[ch].begin(), adj[ch].end(), node) == adj[ch].end()) {
                all_direct = false;
                break;
            }
        if (all_direct) return node;
    }
    return none;
}

/** Port of `buildWorkflow` (with `buildAdjacency` and `addPrecedences`). */
template <class T>
workflow::Workflow<T> build_workflow(const json& data, const std::string& wf_name,
                                     const WfCommonsOptions& opt) {
    typedef workflow::Workflow<T> W;
    W wf(wf_name);
    const json& jwf = data.at("workflow");
    const bool legacy = !has(jwf, "specification");
    const json& tasks = legacy ? jwf.at("tasks") : jwf.at("specification").at("tasks");
    if (!tasks.is_array())
        throw InputError("WfCommonsLoader: the workflow's 'tasks' must be an array of task objects");
    const std::size_t n = tasks.size();

    std::vector<std::string> ids(n);
    for (std::size_t i = 0; i < n; ++i) ids[i] = task_id(tasks[i], i + 1);

    // Execution data, keyed by task id.
    std::map<std::string, const json*> exec_map;
    if (opt.use_execution_data) {
        if (legacy) {
            for (std::size_t i = 0; i < n; ++i) exec_map[ids[i]] = &tasks[i];
        } else if (has(jwf, "execution") && has(jwf.at("execution"), "tasks")) {
            for (const json& et : jwf.at("execution").at("tasks"))
                if (has(et, "id")) exec_map[as_text(et.at("id"))] = &et;
        }
    }

    // Phase 1: activities.
    std::map<std::string, std::size_t> idx_of;
    for (std::size_t i = 0; i < n; ++i) {
        const std::string& id = ids[i];
        double runtime = opt.default_runtime;
        const auto it = exec_map.find(id);
        const json* exec = it == exec_map.end() ? nullptr : it->second;
        if (exec && has(*exec, "runtimeInSeconds") && exec->at("runtimeInSeconds").is_number())
            runtime = exec->at("runtimeInSeconds").get<double>();
        wf.add_activity(id, fit_distribution<T>(runtime, opt));
        idx_of[id] = i;
        if (opt.store_metadata) wf.activity(id).set_metadata(extract_metadata(tasks[i], id, exec));
    }

    // Children named by task name, for the fallback (a name shared by two tasks is ambiguous).
    std::map<std::string, std::size_t> name_idx;
    std::map<std::string, int> name_count;
    if (opt.resolve_children_by_name)
        for (std::size_t i = 0; i < n; ++i)
            if (has(tasks[i], "name") && tasks[i].at("name").is_string()) {
                const std::string nm = tasks[i].at("name").get<std::string>();
                name_idx[nm] = i;
                ++name_count[nm];
            }

    // Phase 2: adjacency, the union of `children` (task -> child) and `parents` (parent -> task),
    // each pair once and the children-derived edges first, as WfCommonsLoader.m builds it.
    std::vector<std::vector<std::size_t>> adj(n);
    std::vector<std::size_t> indeg(n, 0), outdeg(n, 0);
    std::vector<std::vector<bool>> seen(n, std::vector<bool>(n, false));
    std::vector<std::string> unresolved;
    for (int pass = 0; pass < 2; ++pass) {
        const char* field = pass == 0 ? "children" : "parents";
        for (std::size_t i = 0; i < n; ++i) {
            if (!has(tasks[i], field)) continue;
            const json& rj = tasks[i].at(field);
            std::vector<std::string> refs;
            if (rj.is_string()) refs.push_back(rj.get<std::string>());
            else if (rj.is_array())
                for (const json& c : rj) refs.push_back(as_text(c));
            for (const std::string& k : refs) {
                std::size_t j = n;
                const auto a = idx_of.find(k);
                if (a != idx_of.end()) j = a->second;
                else if (opt.resolve_children_by_name) {
                    const auto b = name_idx.find(k);
                    if (b != name_idx.end() && name_count[k] == 1) j = b->second;
                }
                if (j == n) {
                    unresolved.push_back(k);
                    continue;
                }
                const std::size_t from = pass == 0 ? i : j, to = pass == 0 ? j : i;
                if (seen[from][to]) continue;
                seen[from][to] = true;
                adj[from].push_back(to);
                ++outdeg[from];
                ++indeg[to];
            }
        }
    }
    if (!unresolved.empty()) {
        std::vector<std::string> shown;
        for (const std::string& u : unresolved)
            if (std::find(shown.begin(), shown.end(), u) == shown.end()) shown.push_back(u);
        std::cerr << "[LINE] Warning: WfCommonsLoader: " << unresolved.size()
                  << " child/parent reference(s) match no task id"
                  << (opt.resolve_children_by_name ? " or unique task name" : "") << " and were dropped: ";
        for (std::size_t u = 0; u < shown.size() && u < 5; ++u) std::cerr << (u ? ", " : "") << shown[u];
        std::cerr << (shown.size() > 5 ? ", ..." : "") << "\n";
    }

    // Phase 3: fork-join pairs, then the remaining edges in series.
    std::vector<std::vector<bool>> processed(n, std::vector<bool>(n, false));
    for (std::size_t i = 0; i < n; ++i) {
        if (outdeg[i] <= 1) continue;
        const std::vector<std::size_t>& children = adj[i];
        const std::size_t join = find_common_join(children, adj, indeg);
        if (join == n) continue;
        std::vector<std::string> posts;
        for (std::size_t c : children) posts.push_back(ids[c]);
        wf.add_precedence(W::AndFork(ids[i], posts));
        wf.add_precedence(W::AndJoin(posts, ids[join]));
        for (std::size_t c : children) {
            processed[i][c] = true;
            processed[c][join] = true;
        }
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j : adj[i])
            if (!processed[i][j]) wf.add_precedence(W::Serial(ids[i], ids[j]));
    return wf;
}

/** Parse a whole file as JSON, naming the file on failure. */
inline json parse_file(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("WfCommonsLoader: cannot open '" + path + "'");
    try {
        json j;
        in >> j;
        return j;
    } catch (const json::parse_error& e) {
        throw InputError("WfCommonsLoader: malformed JSON in '" + path + "': " + e.what());
    }
}

}  // namespace wfcommons_detail

/**
 * Port of `WfCommonsLoader.loadFromStruct(data, options)`, on a parsed document.
 * The name defaults to `struct_input`, as in the reference.
 */
template <class T>
workflow::Workflow<T> wfcommons_load_from_json(const nlohmann::json& data,
                                               const WfCommonsOptions& options = WfCommonsOptions()) {
    wfcommons_detail::validate_schema(data);
    const std::string nm = wfcommons_detail::extract_name(data, "struct_input");
    return wfcommons_detail::build_workflow<T>(data, nm, options);
}

/** Port of `WfCommonsLoader.load(jsonFile, options)`. */
template <class T>
workflow::Workflow<T> wfcommons_load(const std::string& path,
                                     const WfCommonsOptions& options = WfCommonsOptions()) {
    const nlohmann::json data = wfcommons_detail::parse_file(path);
    wfcommons_detail::validate_schema(data);
    const std::string nm = wfcommons_detail::extract_name(data, path);
    return wfcommons_detail::build_workflow<T>(data, nm, options);
}

namespace wfcommons_detail {

/** The first executable `name` on PATH, empty when there is none. */
inline std::string which_on_path(const std::string& name) {
    const char* env = std::getenv("PATH");
    const std::string path = env ? env : "";
    std::size_t b = 0;
    while (b <= path.size()) {
        const std::size_t e = path.find(':', b);
        const std::string dir = path.substr(b, e == std::string::npos ? std::string::npos : e - b);
        if (!dir.empty() && ::access((dir + "/" + name).c_str(), X_OK) == 0) return dir + "/" + name;
        if (e == std::string::npos) break;
        b = e + 1;
    }
    return std::string();
}

/**
 * GET an https:// URL with curl, else wget, into memory. Redirects are followed but only to https, and an HTTP
 * error status fails the fetch; the tool's own message is quoted in the `InputError`.
 */
inline std::string fetch_https(const std::string& url, int timeout_ms) {
    const int secs = timeout_ms > 0 ? (timeout_ms + 999) / 1000 : 0;
    util::TempDir dir("wfcommons");
    const std::string dest = dir.file("workflow.json");
    std::vector<std::string> argv;
    const std::string curl = which_on_path("curl");
    if (!curl.empty()) {
        argv = {curl, "-sS", "-fL", "--proto", "=https", "--proto-redir", "=https", "--tlsv1.2", "-o", dest};
        if (secs > 0) {
            argv.push_back("--max-time");
            argv.push_back(std::to_string(secs));
        }
    } else {
        const std::string wget = which_on_path("wget");
        if (wget.empty())
            throw UnsupportedError("WfCommonsLoader: " + url + " is https, which this port fetches with curl or "
                                   "wget, and neither is on PATH. Install one, or download the file and read it "
                                   "with wfcommons_load");
        argv = {wget, "-q", "--https-only", "--tries=1", "-O", dest};
        if (secs > 0) argv.push_back("--timeout=" + std::to_string(secs));
    }
    argv.push_back(url);
    const util::ProcResult r = util::capture(argv, secs > 0 ? secs + 5 : 0, true);
    if (r.timedOut) throw InputError("WfCommonsLoader: GET " + url + " did not finish within " +
                                     std::to_string(secs) + " s");
    if (r.exitCode != 0)
        throw InputError("WfCommonsLoader: GET " + url + " failed (" + argv[0] + " exit " +
                         std::to_string(r.exitCode) + "): " + util::trim(r.out));
    std::ifstream in(dest.c_str(), std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace wfcommons_detail

/**
 * Port of `WfCommonsLoader.loadFromUrl(urlString, options)`.
 *
 * http:// goes over `line::http`, where a non-200 answer, including a redirect,
 * is an `InputError` naming the status. https:// goes through curl or wget
 * (`fetch_https`), following redirects that stay on https.
 *
 * @param timeout_ms read timeout for the fetch
 */
template <class T>
workflow::Workflow<T> wfcommons_load_from_url(const std::string& url,
                                              const WfCommonsOptions& options = WfCommonsOptions(),
                                              int timeout_ms = 60000) {
    std::string body;
    if (url.compare(0, 8, "https://") == 0) {
        body = wfcommons_detail::fetch_https(url, timeout_ms);
    } else {
        const http::Response resp = http::get(url, timeout_ms);
        if (resp.status != 200)
            throw InputError("WfCommonsLoader: GET " + url + " answered HTTP " +
                             std::to_string(resp.status) + " (redirects are not followed)");
        body = resp.body;
    }
    nlohmann::json data;
    try {
        data = nlohmann::json::parse(body);
    } catch (const nlohmann::json::parse_error& e) {
        throw InputError("WfCommonsLoader: malformed JSON from " + url + ": " + e.what());
    }
    wfcommons_detail::validate_schema(data);
    const std::string nm = wfcommons_detail::extract_name(data, wfcommons_detail::file_stem(url));
    return wfcommons_detail::build_workflow<T>(data, nm, options);
}

/**
 * Port of `WfCommonsLoader.validateFile(jsonFile)`: true when the file parses and
 * passes the schema checks. Never throws.
 */
inline bool wfcommons_validate(const std::string& path) {
    try {
        wfcommons_detail::validate_schema(wfcommons_detail::parse_file(path));
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_WFCOMMONS_LOADER_H
