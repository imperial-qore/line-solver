/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_WRAPPERS_JMT_JMT_LOGS_H
#define LINE_SOLVERS_WRAPPERS_JMT_JMT_LOGS_H

/**
 * The log-driven half of `SolverJMT`: `linkAndLog`, `parseLogs`,
 * `parseTranState`, `parseTranRespT`, `sampleAggr`, `sampleSysAggr`,
 * `getProbAggr`, `getProbSysAggr`, `getCdfRespT` and `getTranProbAggr`.
 *
 * WHY THERE IS A SECOND MODEL. JMT reports MEANS, not trajectories; everything
 * here recovers a trajectory instead, by rebuilding the model with a Logger on
 * each side of every node of interest, running the simulation, and reading the
 * arrival and departure CSV files back. `jmt_link_and_log` is that rebuild, and
 * it is the port of `@@MNetwork/linkAndLog.m`: a job entering node i now
 * crosses `Arv_i` first and leaves through `Dep_i`, so the two files bracket
 * every passage through i.
 *
 * THE LOGGERS CHANGE THE TOPOLOGY, NOT THE MODEL. They hold no jobs and route
 * with probability one, so the stochastic complement that removes a Router
 * removes them too and the stationary law is unchanged; what they add is the
 * event stream. That is why the sample path this returns is the sample path of
 * the original model and not of an instrumented approximation of it.
 *
 * THE TARGET STATE NOW TRAVELS. `getProbAggr` and `getProbSysAggr` weigh the
 * trajectory against `sn.state{isf}`, the model's CURRENT state, and this
 * header used to refuse both because `qn::NetworkStruct` carried no such thing.
 * It does now: every writer emits the (`stateSpace`, `statePrior`) pair for each
 * stateful node that carries a state -- which is what `setState` and
 * `initFromMarginal` leave -- and `sn_declared_marginal` decodes that pair back
 * into the per-class job counts the two getters compare against, falling back
 * per station to the default marking where nothing was declared.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "line/api/sn/sn_state.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/wrappers/jmt/solver_jmt.h"
#include "line/util/error.h"
#include "line/util/tempdir.h"

namespace line {
namespace jmt {

/** The event kinds a JMT log carries, MATLAB `EventType`. */
enum class JmtEventType { INIT, ARV, DEP };

/** One logged event: when, of which kind, in which class, for which job. */
struct JmtEvent {
    double t = 0.0;
    JmtEventType type = JmtEventType::INIT;
    std::size_t node = 0;   ///< 1-based node index
    std::size_t cls = 0;    ///< 1-based class index, 0 at INIT
    double job = 0.0;       ///< JMT's job id, -1 at INIT
};

/**
 * Port of `@@MNetwork/linkAndLog.m`: the model with an arrival and a departure
 * Logger around every logged node.
 *
 * @param sn             the model to instrument
 * @param is_node_logged (nnodes) which nodes get the Logger pair
 * @param log_path       the directory the CSV files are written to
 *
 * The new node order is the reference's: the original nodes, then the arrival
 * Loggers in node order, then the departure Loggers. A SOURCE AND A SINK ARE
 * NEVER LOGGED -- an arrival Logger before a Source has nothing to observe and
 * a departure Logger after a Sink is unreachable -- and the reference drops
 * them from the request with a warning; this port drops them silently, since
 * the caller's `is_node_logged` is derived from the metric handles and never
 * asks for either on purpose.
 */
template <class T>
qn::NetworkStruct<T> jmt_link_and_log(const qn::NetworkStruct<T>& sn,
                                      const std::vector<bool>& is_node_logged,
                                      const std::string& log_path) {
    const std::size_t N = sn.nodes.size(), K = sn.classes.size();
    if (is_node_logged.size() != N)
        throw InputError("linkAndLog: the isNodeLogged vector does not match the node count");
    std::vector<bool> logged = is_node_logged;
    for (std::size_t i = 0; i < N; ++i)
        if (sn.nodes[i].nodetype == lang::NodeType::Source ||
            sn.nodes[i].nodetype == lang::NodeType::Sink)
            logged[i] = false;

    qn::NetworkStruct<T> out = sn;
    out.log_path = log_path;
    // The routing is rebuilt from scratch below, so every derived table has to
    // go with it; `refresh_struct` at the end rederives them all.
    out.P.clear();
    out.Peff.clear();

    std::vector<std::size_t> arv(N + 1, 0), dep(N + 1, 0);
    for (std::size_t i = 1; i <= N; ++i) {
        if (!logged[i - 1]) continue;
        const std::size_t nd = out.add_node("Arv_" + sn.nodes[i - 1].name, lang::NodeType::Logger,
                                            false);
        // PROB, NOT RAND, AND IT IS THE SAMPLE PATH. A Logger has ONE
        // destination, so the two strategies describe the same routing -- but
        // JSIM spells them differently (an EmpiricalStrategy carrying
        // probability 1, against a RandomStrategy) and draws from the stream
        // differently, so a model whose loggers route RAND walks a DIFFERENT
        // path from the reference's at the same seed. The reference re-links
        // the logged network (`linkAndLog` calls `link`), which sets PROB
        // everywhere, and on cdf_respt_open_twoclasses this one difference was
        // the whole of the C++ row's 5% gap: Queue1 mean response time 2.041
        // against MATLAB's 1.9358 at seed 23000. Same trap as building a twin
        // from a RoutingMatrix where the reference used addLink.
        out.nodes[nd - 1].routing.assign(K, lang::RoutingStrategy::PROB);
        out.nodes[nd - 1].routing_weights.assign(K, std::map<std::size_t, double>());
        out.nodes[nd - 1].routing_param.assign(K, 0);
        out.nodes[nd - 1].logger.file_name = sn.nodes[i - 1].name + "-Arv.csv";
        out.nodes[nd - 1].logger.file_path = log_path;
        arv[i] = nd;
    }
    for (std::size_t i = 1; i <= N; ++i) {
        if (!logged[i - 1]) continue;
        const std::size_t nd = out.add_node("Dep_" + sn.nodes[i - 1].name, lang::NodeType::Logger,
                                            false);
        // PROB for the same reason as the arrival logger above.
        out.nodes[nd - 1].routing.assign(K, lang::RoutingStrategy::PROB);
        out.nodes[nd - 1].routing_weights.assign(K, std::map<std::size_t, double>());
        out.nodes[nd - 1].routing_param.assign(K, 0);
        out.nodes[nd - 1].logger.file_name = sn.nodes[i - 1].name + "-Dep.csv";
        out.nodes[nd - 1].logger.file_path = log_path;
        dep[i] = nd;
    }

    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);
    for (const auto& kv : sn.P) {
        const std::size_t r = kv.first.first, s = kv.first.second;
        const Matrix<T>& B = kv.second;
        for (std::size_t i = 1; i <= N && i <= B.rows(); ++i)
            for (std::size_t j = 1; j <= N && j <= B.cols(); ++j) {
                if (B(i - 1, j - 1) == zero) continue;
                const std::size_t from = logged[i - 1] ? dep[i] : i;
                const std::size_t to = logged[j - 1] ? arv[j] : j;
                out.set_route(r, s, from, to, B(i - 1, j - 1));
            }
    }
    // The two unit arcs that put the loggers on the path: Arv_i -> i -> Dep_i,
    // both class preserving, since a Logger switches nothing.
    for (std::size_t i = 1; i <= N; ++i) {
        if (!logged[i - 1]) continue;
        for (std::size_t r = 1; r <= K; ++r) {
            out.set_route(r, r, arv[i], i, one);
            out.set_route(r, r, i, dep[i], one);
        }
    }
    out.refresh_struct();
    return out;
}

namespace detail {

/** Split a line on `;`, the JSIM header's `logDelimiter`. */
inline std::vector<std::string> split_semi(const std::string& line) {
    std::vector<std::string> out;
    std::size_t b = 0;
    while (true) {
        const std::size_t e = line.find(';', b);
        out.push_back(line.substr(b, e == std::string::npos ? std::string::npos : e - b));
        if (e == std::string::npos) break;
        b = e + 1;
    }
    return out;
}

}  // namespace detail

/** One arrival or departure CSV file, as three parallel columns. */
struct JmtLogFile {
    std::vector<double> ts;
    std::vector<double> job;
    std::vector<std::string> cls;
};

/**
 * Read one JMT log CSV.
 *
 * THE COLUMN LAYOUT IS THE REFERENCE'S, not inferred from the header: the
 * loggers this port installs carry the default flag set of `Logger.m`
 * (timestamp, job id and job class on; wall-clock start time, logger name and
 * the two inter-departure columns off), and `parseLogs.m` reads columns 2, 3
 * and 4 of a six-column line after one header row. Sniffing the header instead
 * would silently accept a file written with other flags and read the wrong
 * column as the timestamp.
 */
inline JmtLogFile jmt_read_log(const std::string& path) {
    std::ifstream f(path.c_str());
    if (!f) throw InputError("SolverJMT: cannot read the JMT log '" + path + "'");
    JmtLogFile out;
    std::string line;
    bool header = true;
    while (std::getline(f, line)) {
        if (header) {
            header = false;
            continue;
        }
        if (line.empty()) continue;
        const std::vector<std::string> col = detail::split_semi(line);
        if (col.size() < 4) continue;
        out.ts.push_back(std::strtod(col[1].c_str(), nullptr));
        out.job.push_back(std::strtod(col[2].c_str(), nullptr));
        out.cls.push_back(col[3]);
    }
    return out;
}

/** The per-class queue-length trajectory of one node, plus its event stream. */
template <class T>
struct JmtNodeTrace {
    std::vector<double> t;                      ///< event times, ascending
    std::vector<std::vector<double>> qlen;      ///< (|t| x nclasses) counts after the event
    std::vector<JmtEvent> event;
};

/**
 * Port of `parseTranState`: the arrival and departure logs merged into a
 * per-class queue-length trajectory.
 *
 * SIMULTANEOUS EVENTS ARE REORDERED, which is the whole reason this is not a
 * merge sort. Two events at the same timestamp involving the SAME job are a
 * pass-through -- the job arrives and departs with no time in between -- and
 * JMT writes them in file order, not in causal order. When the previous event
 * of that job was of the same kind, the reference swaps the event with the next
 * one involving the job, which restores the alternation; without it the
 * cumulative sum below goes negative and the trajectory reports a queue with
 * -1 jobs in it.
 */
template <class T>
JmtNodeTrace<T> jmt_parse_tran_state(const JmtLogFile& arv, const JmtLogFile& dep,
                                     const std::vector<std::size_t>& class_of_arv,
                                     const std::vector<std::size_t>& class_of_dep,
                                     const std::vector<double>& node_preload) {
    const std::size_t K = node_preload.size();
    struct Row {
        double t;
        std::vector<double> delta;
        JmtEventType type;
        std::size_t cls;
        double job;
    };
    std::vector<Row> rows;
    rows.reserve(arv.ts.size() + dep.ts.size());
    for (std::size_t i = 0; i < arv.ts.size(); ++i) {
        Row r;
        r.t = arv.ts[i];
        r.delta.assign(K, 0.0);
        if (class_of_arv[i] >= 1 && class_of_arv[i] <= K) r.delta[class_of_arv[i] - 1] = +1.0;
        r.type = JmtEventType::ARV;
        r.cls = class_of_arv[i];
        r.job = arv.job[i];
        rows.push_back(r);
    }
    for (std::size_t i = 0; i < dep.ts.size(); ++i) {
        Row r;
        r.t = dep.ts[i];
        r.delta.assign(K, 0.0);
        if (class_of_dep[i] >= 1 && class_of_dep[i] <= K) r.delta[class_of_dep[i] - 1] = -1.0;
        r.type = JmtEventType::DEP;
        r.cls = class_of_dep[i];
        r.job = dep.job[i];
        rows.push_back(r);
    }
    // `sortrows` on the timestamp is STABLE in MATLAB, so arrivals keep their
    // file order among themselves and ahead of the departures of the same
    // instant; std::stable_sort reproduces that.
    std::stable_sort(rows.begin(), rows.end(),
                     [](const Row& a, const Row& b) { return a.t < b.t; });

    for (std::size_t ev = 1; ev < rows.size(); ++ev) {
        if (rows[ev].t != rows[ev - 1].t) continue;  // not an instantaneous pair
        const double j = rows[ev].job;
        std::size_t prev = 0;
        bool has_prev = false;
        for (std::size_t k = ev; k > 0; --k)
            if (rows[k - 1].job == j) {
                prev = k - 1;
                has_prev = true;
                break;
            }
        if (!has_prev) continue;
        if (rows[prev].type != rows[ev].type) continue;
        std::size_t next = 0;
        bool has_next = false;
        for (std::size_t k = ev + 1; k < rows.size(); ++k)
            if (rows[k].job == j) {
                next = k;
                has_next = true;
                break;
            }
        if (!has_next) continue;
        std::swap(rows[ev], rows[next]);
    }

    JmtNodeTrace<T> out;
    out.t.push_back(0.0);
    out.qlen.push_back(node_preload);
    JmtEvent init;
    init.t = 0.0;
    init.type = JmtEventType::INIT;
    init.job = -1.0;
    out.event.push_back(init);
    std::vector<double> run = node_preload;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        for (std::size_t r = 0; r < K; ++r) run[r] += rows[i].delta[r];
        out.t.push_back(rows[i].t);
        out.qlen.push_back(run);
        JmtEvent e;
        e.t = rows[i].t;
        e.type = rows[i].type;
        e.cls = rows[i].cls;
        e.job = rows[i].job;
        out.event.push_back(e);
    }
    return out;
}

/**
 * Port of `parseTranRespT`: the per-class response-time samples of one node.
 *
 * A job's passages are recovered by pairing its events IN TIME ORDER: the first
 * arrival, the matching departure, and so on. A job may pass through the same
 * node several times -- a closed model does nothing else -- so the pairing is
 * over the whole per-job sequence and not just its first and last event, and an
 * unmatched trailing arrival (the job was still there when the run ended, or it
 * was dropped) is discarded rather than paired with the run's end.
 */
inline std::map<std::size_t, std::vector<double>> jmt_parse_tran_resp_t(
    const JmtLogFile& arv, const JmtLogFile& dep, const std::vector<std::size_t>& class_of_arv,
    const std::vector<std::size_t>& class_of_dep) {
    struct Ev {
        double t;
        bool is_arv;
        std::size_t cls;
    };
    std::map<double, std::vector<Ev>> by_job;
    for (std::size_t i = 0; i < arv.ts.size(); ++i)
        by_job[arv.job[i]].push_back(Ev{arv.ts[i], true, class_of_arv[i]});
    for (std::size_t i = 0; i < dep.ts.size(); ++i)
        by_job[dep.job[i]].push_back(Ev{dep.ts[i], false, class_of_dep[i]});

    std::map<std::size_t, std::vector<double>> out;
    for (auto& kv : by_job) {
        std::vector<Ev>& evs = kv.second;
        std::stable_sort(evs.begin(), evs.end(),
                         [](const Ev& a, const Ev& b) { return a.t < b.t; });
        std::size_t i = 0;
        while (i + 1 < evs.size()) {
            if (!evs[i].is_arv) {
                ++i;  // a departure with no arrival before it: the preload
                continue;
            }
            if (evs[i + 1].is_arv) {
                ++i;  // two arrivals in a row: the first was dropped
                continue;
            }
            // The CLASS OF THE ARRIVAL owns the sample: a job that switched
            // class while inside the node arrived in the first one, and that is
            // the class whose response time this passage measures.
            out[evs[i].cls].push_back(evs[i + 1].t - evs[i].t);
            i += 2;
        }
    }
    return out;
}

/**
 * Port of `parseLogs`: read every logged node's CSV pair.
 *
 * @param sn        the INSTRUMENTED struct `jmt_link_and_log` returned
 * @param orig      the original struct, whose node names the files are keyed by
 * @param log_path  where the loggers wrote
 *
 * The preload of a node is its `initial_marginal`, the same value the JSIM
 * `preload` block carried, because the logs record only the CROSSINGS: a job
 * that was already at the node when the run started never arrives, so without
 * the preload the cumulative sum starts from zero and every count is short by
 * the initial population.
 */
template <class T>
std::map<std::size_t, JmtNodeTrace<T>> jmt_parse_logs(const qn::NetworkStruct<T>& orig,
                                                      const std::vector<bool>& is_node_logged,
                                                      const std::string& log_path) {
    std::map<std::string, std::size_t> class_of;
    for (std::size_t r = 1; r <= orig.nclasses; ++r) class_of[orig.classes[r - 1].name] = r;

    std::map<std::size_t, JmtNodeTrace<T>> out;
    for (std::size_t ind = 1; ind <= orig.nodes.size(); ++ind) {
        if (!is_node_logged[ind - 1]) continue;
        const std::string base = log_path + "/" + orig.nodes[ind - 1].name;
        const std::string fa = base + "-Arv.csv", fd = base + "-Dep.csv";
        if (!detail::is_file(fa) || !detail::is_file(fd)) continue;
        const JmtLogFile arv = jmt_read_log(fa), dep = jmt_read_log(fd);
        std::vector<std::size_t> ca(arv.cls.size(), 0), cd(dep.cls.size(), 0);
        for (std::size_t i = 0; i < arv.cls.size(); ++i) {
            const auto it = class_of.find(arv.cls[i]);
            ca[i] = it == class_of.end() ? 0 : it->second;
        }
        for (std::size_t i = 0; i < dep.cls.size(); ++i) {
            const auto it = class_of.find(dep.cls[i]);
            cd[i] = it == class_of.end() ? 0 : it->second;
        }
        std::vector<double> preload(orig.nclasses, 0.0);
        const std::size_t ist = orig.nodes[ind - 1].station;
        if (ist != 0) {
            const auto im = orig.initmarking.find(ind);
            if (im != orig.initmarking.end()) {
                for (std::size_t r = 0; r < orig.nclasses && r < im->second.size(); ++r)
                    preload[r] = num_traits<T>::to_double(im->second[r]);
            } else {
                for (std::size_t r = 0; r < orig.nclasses; ++r) {
                    const double n = orig.classes[r].population;
                    if (std::isfinite(n) && orig.classes[r].refstat == ist) preload[r] = n;
                }
            }
        }
        out[ind] = jmt_parse_tran_state<T>(arv, dep, ca, cd, preload);
    }
    return out;
}

/**
 * Port of `sampleAggr`: the queue-length trajectory of one node.
 *
 * The model is instrumented, simulated once, and the node's CSV pair read back.
 * `num_events` bounds how much of the trajectory is returned, as the reference
 * does; JMT CANNOT BE ASKED FOR A NUMBER OF EVENTS AT ONE NODE -- `maxEvents`
 * is global -- so the count is a truncation of what the run produced and not a
 * stopping rule, which is exactly the caveat the reference warns about.
 */
template <class T>
JmtNodeTrace<T> jmt_sample_aggr(const qn::NetworkStruct<T>& sn, std::size_t node,
                                std::size_t num_events, const JmtOptions& opt) {
    if (node == 0 || node > sn.nodes.size())
        throw InputError("SolverJMT: sampleAggr was given a node index out of range");
    std::vector<bool> logged(sn.nodes.size(), false);
    logged[node - 1] = true;
    util::TempDir logs("jmtlog");
    if (opt.keep) logs.keep();
    const qn::NetworkStruct<T> inst = jmt_link_and_log(sn, logged, logs.path());
    solver_jmt_run_analyzer(inst, opt);
    const std::map<std::size_t, JmtNodeTrace<T>> traces = jmt_parse_logs(sn, logged, logs.path());
    const auto it = traces.find(node);
    if (it == traces.end())
        throw NumericError(
            "SolverJMT: the simulation produced no log for node '" + sn.nodes[node - 1].name +
            "'; the run has likely failed before any job crossed it");
    JmtNodeTrace<T> tr = it->second;
    if (num_events > 0 && tr.t.size() > num_events + 1) {
        tr.t.resize(num_events + 1);
        tr.qlen.resize(num_events + 1);
        std::vector<JmtEvent> ev;
        for (std::size_t i = 0; i < tr.event.size(); ++i)
            if (tr.event[i].t <= tr.t.back()) ev.push_back(tr.event[i]);
        tr.event = ev;
    }
    return tr;
}

/** The system trajectory: one per-class block per station, on a common grid. */
template <class T>
struct JmtSysTrace {
    std::vector<double> t;
    /** state[ist-1] is (|t| x nclasses); a Source's block is empty. */
    std::vector<std::vector<std::vector<double>>> state;
    std::vector<JmtEvent> event;
};

/**
 * Port of `sampleSysAggr`: every station's trajectory on one time grid.
 *
 * The per-station traces have their OWN event times, so they are resampled onto
 * the union grid with a PREVIOUS-value hold -- a queue length is constant
 * between its own events, so holding is exact rather than an interpolation. The
 * grid is truncated at the earliest station's last event: past that point one
 * station has no data, and extending its last value would report a queue that
 * stopped changing rather than one that stopped being observed.
 */
template <class T>
JmtSysTrace<T> jmt_sample_sys_aggr(const qn::NetworkStruct<T>& sn, std::size_t num_events,
                                   const JmtOptions& opt) {
    std::vector<bool> logged(sn.nodes.size(), false);
    for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
        const std::size_t ind = sn.station_to_node[ist - 1];
        if (sn.nodes[ind - 1].nodetype != lang::NodeType::Source) logged[ind - 1] = true;
    }
    util::TempDir logs("jmtlog");
    if (opt.keep) logs.keep();
    const qn::NetworkStruct<T> inst = jmt_link_and_log(sn, logged, logs.path());
    JmtOptions o = opt;
    o.method = "jsim";  // never `replication`, which is composed FROM this
    solver_jmt_run_analyzer(inst, o);
    const std::map<std::size_t, JmtNodeTrace<T>> traces = jmt_parse_logs(sn, logged, logs.path());

    JmtSysTrace<T> out;
    std::set<double> grid;
    double tmax = std::numeric_limits<double>::infinity();
    for (const auto& kv : traces) {
        if (kv.second.t.empty()) continue;
        tmax = std::min(tmax, kv.second.t.back());
        for (double t : kv.second.t) grid.insert(t);
    }
    for (double t : grid)
        if (t <= tmax) out.t.push_back(t);
    if (num_events > 0 && out.t.size() > num_events) out.t.resize(num_events);

    out.state.assign(sn.nstations, std::vector<std::vector<double>>());
    for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
        const std::size_t ind = sn.station_to_node[ist - 1];
        const auto it = traces.find(ind);
        if (it == traces.end()) continue;
        const JmtNodeTrace<T>& tr = it->second;
        std::vector<std::vector<double>> block(out.t.size(),
                                               std::vector<double>(sn.nclasses, 0.0));
        std::size_t k = 0;
        for (std::size_t g = 0; g < out.t.size(); ++g) {
            while (k + 1 < tr.t.size() && tr.t[k + 1] <= out.t[g]) ++k;
            if (k < tr.qlen.size()) block[g] = tr.qlen[k];
        }
        out.state[ist - 1] = block;
    }
    for (const auto& kv : traces)
        for (std::size_t i = 0; i < kv.second.event.size(); ++i) {
            JmtEvent e = kv.second.event[i];
            e.node = kv.first;
            if (e.t <= (out.t.empty() ? 0.0 : out.t.back())) out.event.push_back(e);
        }
    std::stable_sort(out.event.begin(), out.event.end(),
                     [](const JmtEvent& a, const JmtEvent& b) { return a.t < b.t; });
    return out;
}

/**
 * Port of `getCdfRespT`: the empirical response-time distribution per
 * (station, class), as the (F, X) pairs `ecdf` returns.
 *
 * Every station with a service process is logged, the model is simulated, and
 * each node's passages are turned into samples. With `seed_from_steady` (the
 * getCdfRespT contract) the model is FIRST solved for its steady-state queue
 * lengths and the logged run starts preloaded at their rounded values, the
 * closed-class remainder on the bottleneck -- the reference's two-run pipeline,
 * which shortens the warmup and is what makes the seeded curve comparable to
 * MATLAB's for the same seed. Without it (the getTranCdfRespT contract) the
 * logged run starts from the model's default initial state, so the samples
 * cover the transient.
 */
template <class T>
std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<double, double>>>
jmt_get_cdf_resp_t(const qn::NetworkStruct<T>& sn, const JmtOptions& opt,
                   bool seed_from_steady = true) {
    std::vector<bool> logged(sn.nodes.size(), false);
    std::vector<bool> cacheclass = io::jmt_cache_classes(sn);
    for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
        const std::size_t ind = sn.station_to_node[ist - 1];
        for (std::size_t r = 1; r <= sn.nclasses; ++r)
            if (io::jmt_metric_enabled(sn, io::JmtMetricKind::RespT, ist, r, cacheclass))
                logged[ind - 1] = true;
    }
    util::TempDir logs("jmtcdf");
    if (opt.keep) logs.keep();
    qn::NetworkStruct<T> inst = jmt_link_and_log(sn, logged, logs.path());
    if (seed_from_steady) {
        // The reference's first run: steady-state queue lengths, floored per
        // (station, class), a closed class's remainder pushed onto its fullest
        // station so the population balances. The rounded counts ride into the
        // logged copy through `initmarking`, the same channel a Place's tokens
        // take, which is what the writer's preload block reads first.
        const JmtResult<T> pre = solver_jmt_run_analyzer(sn, opt);
        std::vector<std::vector<double>> n(sn.nstations,
                                           std::vector<double>(sn.nclasses, 0.0));
        for (std::size_t r = 1; r <= sn.nclasses; ++r) {
            double tot = 0.0, vmax = -1.0;
            std::size_t imax = 1;
            for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
                const double q = std::floor(std::max(
                    0.0, num_traits<T>::to_double(pre.avg.QN(ist - 1, r - 1))));
                n[ist - 1][r - 1] = q;
                tot += q;
                if (q > vmax) {
                    vmax = q;
                    imax = ist;
                }
            }
            const double njobs = sn.classes[r - 1].population;
            if (std::isfinite(njobs) && tot < njobs)
                n[imax - 1][r - 1] += njobs - tot;
        }
        for (std::size_t ist = 1; ist <= inst.nstations && ist <= sn.nstations; ++ist) {
            const std::size_t ind = inst.station_to_node[ist - 1];
            const lang::NodeType ty = inst.nodes[ind - 1].nodetype;
            if (ty == lang::NodeType::Source || ty == lang::NodeType::Join) continue;
            std::vector<T> row(sn.nclasses);
            for (std::size_t r = 0; r < sn.nclasses; ++r)
                row[r] = num_traits<T>::from_double(n[ist - 1][r]);
            inst.initmarking[ind] = row;
        }
    }
    solver_jmt_run_analyzer(inst, opt);

    std::map<std::string, std::size_t> class_of;
    for (std::size_t r = 1; r <= sn.nclasses; ++r) class_of[sn.classes[r - 1].name] = r;

    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<double, double>>> out;
    for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
        const std::size_t ind = sn.station_to_node[ist - 1];
        if (!logged[ind - 1]) continue;
        const std::string base = logs.path() + "/" + sn.nodes[ind - 1].name;
        if (!detail::is_file(base + "-Arv.csv") || !detail::is_file(base + "-Dep.csv")) continue;
        const JmtLogFile arv = jmt_read_log(base + "-Arv.csv");
        const JmtLogFile dep = jmt_read_log(base + "-Dep.csv");
        std::vector<std::size_t> ca(arv.cls.size(), 0), cd(dep.cls.size(), 0);
        for (std::size_t i = 0; i < arv.cls.size(); ++i) {
            const auto it = class_of.find(arv.cls[i]);
            ca[i] = it == class_of.end() ? 0 : it->second;
        }
        for (std::size_t i = 0; i < dep.cls.size(); ++i) {
            const auto it = class_of.find(dep.cls[i]);
            cd[i] = it == class_of.end() ? 0 : it->second;
        }
        const std::map<std::size_t, std::vector<double>> samples =
            jmt_parse_tran_resp_t(arv, dep, ca, cd);
        for (const auto& kv : samples) {
            if (kv.first == 0 || kv.second.empty()) continue;
            std::vector<double> x = kv.second;
            std::sort(x.begin(), x.end());
            // `ecdf`'s convention: the first row is (0, x(1)) and the step at
            // each distinct value carries the mass of its ties.
            std::vector<std::pair<double, double>> fx;
            fx.push_back(std::make_pair(0.0, x.front()));
            std::size_t i = 0;
            while (i < x.size()) {
                std::size_t j = i;
                while (j + 1 < x.size() && x[j + 1] == x[i]) ++j;
                fx.push_back(std::make_pair(static_cast<double>(j + 1) /
                                                static_cast<double>(x.size()),
                                            x[i]));
                i = j + 1;
            }
            out[std::make_pair(ist, kv.first)] = fx;
        }
    }
    return out;
}

/**
 * Port of `getTranProbAggr`: the transient distribution of one station's
 * aggregate state, estimated over `replications` independent runs.
 *
 * Each run is a fresh seed, `seed + it`, and the runs are resampled onto the
 * union of their event times before being counted, so a state's probability at
 * time t is the fraction of replications occupying it at t. A FINITE HORIZON IS
 * REQUIRED: with an unbounded one every replication ends at a different time
 * and the union grid is not a grid the estimate is defined on -- which is
 * exactly the error the reference raises.
 */
template <class T>
std::pair<std::vector<double>, std::vector<std::vector<double>>> jmt_get_tran_prob_aggr(
    const qn::NetworkStruct<T>& sn, std::size_t station, std::size_t replications,
    const JmtOptions& opt, std::vector<std::vector<double>>& states_out) {
    if (!std::isfinite(opt.max_simulated_time))
        throw InputError(
            "SolverJMT: getTranProbAggr requires a finite time span, e.g. "
            "options.timespan = [0, T]");
    if (station == 0 || station > sn.nstations)
        throw InputError("SolverJMT: getTranProbAggr was given a station index out of range");
    const std::size_t ind = sn.station_to_node[station - 1];
    if (sn.nodes[ind - 1].nodetype == lang::NodeType::Source)
        throw InputError("SolverJMT: getTranProbAggr does not apply to a Source");

    std::vector<JmtSysTrace<T>> runs;
    std::set<double> grid;
    for (std::size_t it = 0; it < replications; ++it) {
        JmtOptions o = opt;
        o.seed = opt.seed + static_cast<long>(it);
        JmtSysTrace<T> tr = jmt_sample_sys_aggr(sn, 0, o);
        for (double t : tr.t) grid.insert(t);
        runs.push_back(tr);
    }
    std::vector<double> tu(grid.begin(), grid.end());

    // The distinct aggregate states seen, in first-appearance order after the
    // resampling; a replication that produced no data contributes none.
    std::map<std::vector<double>, std::size_t> index;
    std::vector<std::vector<double>> states;
    std::vector<std::vector<std::size_t>> occupied(runs.size(),
                                                   std::vector<std::size_t>(tu.size(), 0));
    for (std::size_t r = 0; r < runs.size(); ++r) {
        const std::vector<std::vector<double>>& block = runs[r].state[station - 1];
        if (block.empty()) continue;
        std::size_t k = 0;
        for (std::size_t g = 0; g < tu.size(); ++g) {
            while (k + 1 < runs[r].t.size() && runs[r].t[k + 1] <= tu[g]) ++k;
            if (k >= block.size()) continue;
            const std::vector<double>& s = block[k];
            auto it = index.find(s);
            if (it == index.end()) {
                index[s] = states.size();
                states.push_back(s);
                it = index.find(s);
            }
            occupied[r][g] = it->second + 1;  // 1-based; 0 marks "no data"
        }
    }
    std::vector<std::vector<double>> pi(tu.size(), std::vector<double>(states.size(), 0.0));
    for (std::size_t r = 0; r < runs.size(); ++r)
        for (std::size_t g = 0; g < tu.size(); ++g)
            if (occupied[r][g] != 0)
                pi[g][occupied[r][g] - 1] += 1.0 / static_cast<double>(runs.size());
    states_out = states;
    return std::make_pair(tu, pi);
}

/**
 * The dwell-time weights of a trajectory sampled at `t`.
 *
 * The reference writes `dt = [diff(t); 0]`: a state observed at t_k is held
 * until t_{k+1}, and the LAST sample carries no weight because the run ends
 * there and how long it would have been held is not observed. Dropping the
 * trailing zero instead -- weighting the last sample by the previous gap -- adds
 * a dwell that was never measured, which is why it is kept as a zero.
 */
inline std::vector<double> jmt_dwell_weights(const std::vector<double>& t) {
    std::vector<double> dt(t.size(), 0.0);
    for (std::size_t i = 0; i + 1 < t.size(); ++i) dt[i] = t[i + 1] - t[i];
    return dt;
}

/** What `jmt_prob_aggr` reports: the system probability and the per-station ones. */
struct JmtProbAggr {
    double sys = 0.0;                 ///< P(the whole network is in the declared state)
    std::vector<double> station;      ///< P(station i holds its declared per-class counts)
    bool sys_seen = false;            ///< whether the joint state occurred at all
    std::vector<bool> station_seen;   ///< the same, per station
};

/**
 * Port of `getProbAggr` and `getProbSysAggr`, both off ONE instrumented run.
 *
 * THEY ARE TIME AVERAGES, NOT SAMPLE FRACTIONS. The trajectory is event-driven,
 * so its samples are not equally spaced; counting them would weigh a state the
 * network leaves at once as heavily as one it sits in. Each sample is weighed by
 * the interval it is held for, which is what makes the estimate converge to the
 * stationary probability.
 *
 * ONE RUN, NOT ONE PER QUESTION. The reference instruments one node at a time
 * (`sampleAggr` logs the node it was asked about) and therefore simulates again
 * for every station; `jmt_sample_sys_aggr` already logs every station on a
 * common grid, so the station marginals are read off the same sample path as the
 * joint probability. That is a DIFFERENT estimator from the reference's -- same
 * quantity, one sample path instead of M independent ones -- and it is the
 * cheaper and the more consistent of the two, since the marginals it reports
 * cannot contradict the joint they were taken from.
 *
 * `target` OVERRIDES THE DECLARED STATE of one station where it is given, which
 * is `getProbAggr`'s second argument; every other station keeps the row
 * `setState`/`initFromMarginal` left, exactly as the reference substitutes into
 * `sn.state{isf}` and leaves the rest alone.
 *
 * A STATE NEVER SEEN IS 0 AND IS REPORTED AS SUCH, with `seen` false beside it.
 * The reference warns and returns zero rather than erroring, because on a
 * simulation a state of small probability legitimately fails to appear in a
 * finite run: that is a statement about the run length, and the caller is the
 * one who can lengthen it.
 */
template <class T>
JmtProbAggr jmt_prob_aggr(const qn::NetworkStruct<T>& sn, const JmtOptions& opt,
                          std::size_t target_station = 0,
                          const std::vector<double>& target = std::vector<double>()) {
    Matrix<T> nirm = api::sn_declared_marginal(sn);
    if (!target.empty()) {
        if (target_station == 0 || target_station > sn.nstations)
            throw InputError("SolverJMT: getProbAggr was given a station index out of range");
        if (target.size() != sn.nclasses)
            throw InputError("SolverJMT: getProbAggr takes one job count per class (" +
                             std::to_string(sn.nclasses) + "), got " +
                             std::to_string(target.size()));
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            nirm(target_station - 1, r) = num_traits<T>::from_double(target[r]);
    }
    std::vector<std::vector<double>> nir(sn.nstations, std::vector<double>(sn.nclasses, 0.0));
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            nir[i][r] = num_traits<T>::to_double(nirm(i, r));

    const JmtSysTrace<T> tr = jmt_sample_sys_aggr(sn, 0, opt);
    const std::vector<double> dt = jmt_dwell_weights(tr.t);

    JmtProbAggr out;
    out.station.assign(sn.nstations, 0.0);
    out.station_seen.assign(sn.nstations, false);
    std::vector<double> hit(sn.nstations, 0.0);
    double sys_hit = 0.0, total = 0.0;
    for (std::size_t g = 0; g < dt.size(); ++g) {
        total += dt[g];
        bool all = true;
        for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
            const std::size_t ind = sn.station_to_node[ist - 1];
            // A Source holds no jobs and the reference's own block for it is not
            // a state; it takes no part in either comparison.
            if (sn.nodes[ind - 1].nodetype == lang::NodeType::Source) continue;
            const std::vector<std::vector<double>>& block = tr.state[ist - 1];
            bool same = g < block.size();
            for (std::size_t r = 0; r < sn.nclasses && same; ++r)
                same = block[g][r] == nir[ist - 1][r];
            if (same) {
                hit[ist - 1] += dt[g];
                if (dt[g] > 0.0) out.station_seen[ist - 1] = true;
            } else {
                all = false;
            }
        }
        if (all) {
            sys_hit += dt[g];
            if (dt[g] > 0.0) out.sys_seen = true;
        }
    }
    // A run that observed no dwell at all has measured nothing, and 0 would read
    // as "the state was never visited" -- a claim the run cannot support.
    if (!(total > 0.0))
        throw NumericError(
            "SolverJMT: the simulation produced no observed dwell time, so no state probability "
            "can be estimated from it");
    out.sys = sys_hit / total;
    for (std::size_t i = 0; i < sn.nstations; ++i) out.station[i] = hit[i] / total;
    return out;
}

/** Transient averages over independent replications, on one time grid. */
template <class T>
struct JmtReplication {
    std::vector<double> t;
    /** QNt[ist-1][r] over `t`. */
    std::vector<std::vector<std::vector<double>>> QNt;
    std::vector<std::vector<std::vector<double>>> UNt;
    std::vector<std::vector<std::vector<double>>> TNt;
    /** Replications that produced a usable trajectory. */
    std::size_t valid = 0;
};

/**
 * Port of the `replication` method of `@@SolverJMT/runAnalyzer.m`.
 *
 * A SINGLE SAMPLE PATH IS NOT THE TRANSIENT MEAN. There is no time-ergodicity at
 * a fixed t, so E[N](t) is estimated by averaging `iter_max` INDEPENDENT
 * replications, each seeded `seed + it`, rather than by reading one trajectory.
 *
 * The replications have their own event grids, so they are resampled onto the
 * union with a previous-value hold and the grid is truncated at the MINIMUM of
 * their last events: past that point one replication has no data, and extending
 * its last value would report a queue that stopped changing rather than one that
 * stopped being observed -- the same rule `jmt_sample_sys_aggr` applies across
 * stations, applied here across seeds.
 *
 * Utilization is `min(n, c)/c` at a finite server and the raw queue length at a
 * delay; throughput follows it as `U*c*mu` and `U*mu`, which is how the
 * reference reads departures off the occupancy rather than differencing the
 * trajectory.
 */
template <class T>
JmtReplication<T> jmt_replication(const qn::NetworkStruct<T>& sn, const JmtOptions& opt) {
    // The predicate `auto_family_refusal` asks, so the gate that decides whether
    // to OFFER 'replication' and this run cannot drift apart.
    {
        const std::string refusal = jmt_method_refusal(sn, "replication", opt);
        if (!refusal.empty()) throw InputError(refusal);
    }
    const int reps = opt.iter_max > 0 ? opt.iter_max : 10;

    std::vector<JmtSysTrace<T>> paths;
    std::set<double> grid;
    double tumax = std::numeric_limits<double>::infinity();
    for (int it = 0; it < reps; ++it) {
        JmtOptions o = opt;
        o.method = "jsim";
        o.seed = opt.seed + it;
        JmtSysTrace<T> tr;
        try {
            tr = jmt_sample_sys_aggr(sn, 0, o);
        } catch (const std::exception&) {
            continue;  // a replication that produced no log contributes nothing
        }
        if (tr.t.empty()) continue;
        tumax = std::min(tumax, tr.t.back());
        for (double tv : tr.t) grid.insert(tv);
        paths.push_back(tr);
    }
    if (paths.empty())
        throw UnsupportedError(
            "SolverJMT: no valid replications produced; the transient averages cannot be computed");

    JmtReplication<T> out;
    out.valid = paths.size();
    for (double tv : grid)
        if (tv <= tumax) out.t.push_back(tv);
    const std::size_t nt = out.t.size();
    const std::size_t M = sn.nstations, K = sn.nclasses;
    out.QNt.assign(M, std::vector<std::vector<double>>(K, std::vector<double>(nt, 0.0)));
    out.UNt.assign(M, std::vector<std::vector<double>>(K, std::vector<double>(nt, 0.0)));
    out.TNt.assign(M, std::vector<std::vector<double>>(K, std::vector<double>(nt, 0.0)));

    const double inv = 1.0 / static_cast<double>(paths.size());
    for (std::size_t ist = 0; ist < M; ++ist) {
        const double c = num_traits<T>::to_double(sn.stations[ist].nservers);
        const bool finite_c = std::isfinite(c) && c > 0.0;
        for (const JmtSysTrace<T>& tr : paths) {
            if (ist >= tr.state.size() || tr.state[ist].empty()) continue;
            const std::vector<std::vector<double>>& block = tr.state[ist];
            std::size_t k = 0;
            for (std::size_t g = 0; g < nt; ++g) {
                while (k + 1 < tr.t.size() && tr.t[k + 1] <= out.t[g]) ++k;
                if (out.t[g] < tr.t.front() || k >= block.size()) continue;
                for (std::size_t r = 0; r < K && r < block[k].size(); ++r) {
                    const double n = block[k][r];
                    out.QNt[ist][r][g] += inv * n;
                    out.UNt[ist][r][g] += inv * (finite_c ? std::min(n, c) / c : n);
                }
            }
        }
        for (std::size_t r = 0; r < K; ++r) {
            const double mu = num_traits<T>::to_double(sn.rates(ist, r));
            const double scale = std::isfinite(mu) ? (finite_c ? c * mu : mu) : 0.0;
            for (std::size_t g = 0; g < nt; ++g) out.TNt[ist][r][g] = out.UNt[ist][r][g] * scale;
        }
    }
    return out;
}

}  // namespace jmt
}  // namespace line

#endif  // LINE_SOLVERS_WRAPPERS_JMT_JMT_LOGS_H
