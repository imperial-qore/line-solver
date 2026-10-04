/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_CODE_GEN_H
#define LINE_IO_CODE_GEN_H

/**
 * @file
 * @ingroup line_io
 * Port of the model -> source code generators of `matlab/src/io/`:
 * `QN2MATLAB.m`, `QN2JAVA.m`, `LQN2JAVA.m`, the generator form of `LQN2MATLAB.m`,
 * and the `LINE2MATLAB.m` / `LINE2JAVA.m` dispatchers. LQN2MATLAB, unlike the
 * others, re-declares the model rather than its struct; see lqn2matlab below.
 *
 * WHAT IS WRITTEN IS THE STRUCT, NOT THE DECLARATIONS. Exactly as in MATLAB the
 * queueing generators read `sn` after refresh: every process is re-expressed
 * from the moments of its Markovian representation (`sn.proc`), so an SCV of
 * at least 0.5 is written as `Exp.fitMean` (SCV 1) or `APH.fitMeanAndSCV`, a
 * smaller one as an Erlang with round(1/SCV) phases, and a process with no
 * representation as `Disabled`. The routing is the node-level `sn.rtnodes`,
 * with a ClassSwitch written as a Router because the switching is already in
 * the matrix. The layered generator reads the `LayeredNetworkStruct`, plus the
 * declared reply activities and AND-join quorums the struct does not keep.
 *
 * NUMERIC TEXT FOLLOWS MATLAB'S fprintf. `%d` of a non-integer double prints in
 * `%e` form there (a routing share of 1/4 is `2.500000e-01`), infinities print
 * as `Inf` and NaN as `NaN`; the helpers below reproduce that, so a line of
 * output is byte for byte the line MATLAB writes for the same struct.
 *
 * WHERE THE C++ OUTPUT DELIBERATELY DIFFERS FROM MATLAB, it is because the
 * MATLAB text is not valid source for its target:
 *   - QN2JAVA: `new Erlang(rate, n)` with an integer phase count and
 *     `new Replayer("f")`, where MATLAB writes `Erlang(%f,%f)` and
 *     `Replayer("f")`, neither of which compiles.
 * LQN2JAVA follows the MATLAB text exactly, distributions included (`javaDist`).
 * A node type the MATLAB generators have no statement for (Cache, Logger,
 * Place, Transition) is refused by name instead of being omitted, since
 * omitting it writes a script that builds a different model.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/io/lqn_json_reader.h"
#include "line/io/network_reader.h"
#include "line/lang/distribution.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/lqn/lqn_writer.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace code_gen_detail {

/** `%f`, `%e` or `%g` as MATLAB's fprintf prints it: C's text, with `Inf` / `NaN` spelt MATLAB's way. */
inline std::string mfmt(const char* spec, double x) {
    if (std::isnan(x)) return "NaN";
    if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
    char buf[64];
    std::snprintf(buf, sizeof(buf), spec, x);
    return std::string(buf);
}

inline std::string fmt_f(double x) { return mfmt("%f", x); }
inline std::string fmt_g(double x) { return mfmt("%g", x); }

/** MATLAB `%d` of a double: the integer when it is one, else MATLAB's `%e` override. */
inline std::string fmt_d(double x) {
    if (std::isnan(x) || std::isinf(x)) return mfmt("%f", x);
    if (x == std::floor(x) && std::fabs(x) < 9.007199254740992e15) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.0f", x == 0.0 ? 0.0 : x);
        return std::string(buf);
    }
    return mfmt("%e", x);
}

inline std::string fmt_d(std::size_t x) { return std::to_string(x); }
inline std::string fmt_d(int x) { return std::to_string(x); }

/** `SchedStrategy.toProperty(SchedStrategy.toText(s))`: the enum constant name. */
inline std::string sched_property(lang::SchedStrategy s) {
    const std::string t = lang::sched_to_text(s);
    if (t == "none")
        throw UnsupportedError("code_gen: a station carries no scheduling strategy to write");
    std::string out;
    for (std::size_t i = 0; i < t.size(); ++i)
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(t[i]))));
    return out;
}

/** `strrep(SchedStrategy.toFeature(s),'_','.')`, the layered generator's spelling. */
inline std::string sched_feature(lang::SchedStrategy s) { return "SchedStrategy." + sched_property(s); }

/** One (station, class) process, reduced to what QN2MATLAB / QN2JAVA decide on. */
struct ProcSpec {
    enum Kind { SKIP, REPLAYER, IMMEDIATE, EXP, APH, ERLANG, DISABLED } kind = SKIP;
    bool arrival = false;  ///< the station is EXT scheduled: setArrival, not setService
    double mean = 0.0;
    double scv = 0.0;
    double nphases = 0.0;
    std::string file;
};

/**
 * The branch of QN2MATLAB's process block for station `i`, class `k` (0-based).
 *
 * The Replayer test is MATLAB's `isprop(station,'serviceProcess')`, a property
 * only Queue and its Delay subclass declare, so a trace at a Source goes down
 * the moment path and is written as the fit of its moments, as in MATLAB.
 */
template <class T>
ProcSpec proc_spec(const qn::NetworkStruct<T>& sn, std::size_t i, std::size_t k) {
    ProcSpec ps;
    const std::size_t nd = sn.station_to_node[i];
    const lang::NodeType nt = sn.nodes[nd - 1].nodetype;
    if (nt == lang::NodeType::Join) return ps;
    ps.arrival = sn.stations[i].sched == lang::SchedStrategy::EXT;
    const lang::Distrib<T>& d = sn.service[i][k];
    const bool has_service_prop = nt == lang::NodeType::Queue || nt == lang::NodeType::Delay;
    if (has_service_prop && !d.disabled && d.type == lang::ProcessType::REPLAYER) {
        if (d.trace_file.empty())
            throw UnsupportedError("code_gen: the Replayer of class '" + sn.classes[k].name + "' at '" +
                                   sn.nodes[nd - 1].name +
                                   "' was built from in-memory samples, so there is no trace file to name");
        ps.kind = ProcSpec::REPLAYER;
        ps.file = d.trace_file;
        return ps;
    }
    if (d.disabled) {  // sn.proc is NaN: map_scv is NaN, so MATLAB reaches the isnan arm
        ps.kind = ProcSpec::DISABLED;
        return ps;
    }
    const mam::Map<T> m = lang::dist_to_map(d);
    ps.scv = num_traits<T>::to_double(mam::map_scv(m));
    ps.mean = num_traits<T>::to_double(mam::map_mean(m));
    if (ps.scv >= 0.5) {
        if (ps.scv == 1.0)
            ps.kind = ps.mean < lang::GlobalConstants::CoarseTol ? ProcSpec::IMMEDIATE : ProcSpec::EXP;
        else
            ps.kind = ProcSpec::APH;
    } else {
        ps.kind = ProcSpec::ERLANG;
        ps.nphases = std::max(1.0, std::round(1.0 / ps.scv));
    }
    return ps;
}

/** `find(sn.fj(:,j))`: the 1-based Fork node the 1-based Join node `j` closes. */
template <class T>
std::size_t fork_of(const qn::NetworkStruct<T>& sn, std::size_t j) {
    for (std::size_t f = 0; f < sn.fj.size(); ++f)
        if (sn.fj[f].second == j) return sn.fj[f].first;
    throw InputError("code_gen: Join '" + sn.nodes[j - 1].name +
                     "' closes no Fork, so the model cannot be written as source");
}

/**
 * `zeroPopRefNode(sn, k)` of QN2MATLAB / QN2JAVA: the NODE index (1-based) of the
 * reference station of zero-population closed class `k` (0-based). It is the
 * refstat of the first populated closed class of k's chain (k's own refstat if
 * there is none), else the first station whose `sn.proc{i}{k}{1}` has a nonzero
 * entry, NaN counting as nonzero under `nnz`, so a disabled pair qualifies.
 */
template <class T>
std::size_t empty_class_ref(const qn::NetworkStruct<T>& sn, std::size_t k) {
    std::size_t ist = 0;
    for (std::size_t c = 0; c < sn.chains.size(); ++c) {
        if (k >= sn.chains[c].size() || !sn.chains[c][k]) continue;
        std::size_t cand = 0;
        for (std::size_t r = 0; r < sn.classes.size() && r < sn.chains[c].size(); ++r) {
            const double n = sn.classes[r].population;
            if (sn.chains[c][r] && n > 0 && std::isfinite(n)) {
                cand = sn.classes[r].refstat;
                break;
            }
        }
        if (cand == 0) cand = sn.classes[k].refstat;
        if (cand >= 1 && cand <= sn.stations.size()) ist = cand;
        break;
    }
    for (std::size_t i = 0; ist == 0 && i < sn.stations.size(); ++i) {
        const lang::Distrib<T>& d = sn.service[i][k];
        if (d.disabled || sn.nodes[sn.station_to_node[i] - 1].nodetype == lang::NodeType::Join) {
            ist = i + 1;
            break;
        }
        const mam::Map<T> m = lang::dist_to_map(d);
        for (std::size_t a = 0; a < m.D0.rows() && ist == 0; ++a)
            for (std::size_t b = 0; b < m.D0.cols(); ++b)
                if (num_traits<T>::to_double(m.D0(a, b)) != 0.0) {
                    ist = i + 1;
                    break;
                }
    }
    return ist > 0 ? sn.station_to_node[ist - 1] : 0;
}

/** Refuses a node type neither MATLAB generator has a statement for. */
template <class T>
void check_node_types(const qn::NetworkStruct<T>& sn, const char* who) {
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        switch (sn.nodes[i].nodetype) {
            case lang::NodeType::Source:
            case lang::NodeType::Delay:
            case lang::NodeType::Queue:
            case lang::NodeType::Router:
            case lang::NodeType::Fork:
            case lang::NodeType::Join:
            case lang::NodeType::Sink:
            case lang::NodeType::ClassSwitch:
                break;
            default:
                throw UnsupportedError(std::string(who) + ": node '" + sn.nodes[i].name + "' is a " +
                                       lang::node_type_to_text(sn.nodes[i].nodetype) +
                                       ", which the generator has no statement for");
        }
    }
}

/** The (k, c, i, m, p) of every positive `sn.rtnodes` entry, in QN2MATLAB's loop order, 0-based.
 * Neither generator writes a Sink row or a closed-class Source row. */
struct Route {
    std::size_t k, c, i, m;
    double p;
};

template <class T>
std::vector<Route> routes(const qn::NetworkStruct<T>& sn) {
    std::vector<Route> out;
    const std::size_t K = sn.classes.size(), I = sn.nodes.size();
    if (sn.rtnodes.rows() != I * K) return out;
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t c = 0; c < K; ++c)
            for (std::size_t i = 0; i < I; ++i)
                for (std::size_t m = 0; m < I; ++m) {
                    const double p = num_traits<T>::to_double(sn.rtnodes(i * K + k, m * K + c));
                    const lang::NodeType nt = sn.nodes[i].nodetype;
                    if (!(p > 0) || nt == lang::NodeType::Sink) continue;
                    if (nt == lang::NodeType::Source && std::isfinite(sn.classes[k].population)) continue;
                    out.push_back(Route{k, c, i, m, p});
                }
    return out;
}

/** Opens `path` for writing or throws, naming it. */
inline std::ofstream open_out(const std::string& path) {
    std::ofstream f(path.c_str());
    if (!f) throw InputError("code_gen: cannot open " + path + " for writing");
    return f;
}

}  // namespace code_gen_detail

// ---------------------------------------------------------------------------
// QN2MATLAB
// ---------------------------------------------------------------------------

/**
 * Port of `QN2MATLAB(model, modelName, fid)`: a MATLAB script that rebuilds the
 * network from its refreshed struct.
 */
template <class T>
void qn2matlab(const qn::NetworkStruct<T>& sn, const std::string& model_name, std::ostream& os) {
    using namespace code_gen_detail;
    using lang::NodeType;
    check_node_types(sn, "QN2MATLAB");
    const std::size_t K = sn.classes.size();
    os << "model = Network('" << model_name << "');\n";
    os << "\n%% Block 1: nodes";
    os << "\n";
    for (std::size_t n = 0; n < sn.nodes.size(); ++n) {
        const std::size_t i = n + 1;
        const qn::NodeDef& nd = sn.nodes[n];
        switch (nd.nodetype) {
            case NodeType::Source:
                os << "node{" << i << "} = Source(model, '" << nd.name << "');\n";
                break;
            case NodeType::Delay:
                os << "node{" << i << "} = DelayStation(model, '" << nd.name << "');\n";
                break;
            case NodeType::Queue: {
                const qn::Station<T>& st = sn.stations[nd.station - 1];
                os << "node{" << i << "} = Queue(model, '" << nd.name << "', SchedStrategy."
                   << sched_property(st.sched) << ");\n";
                if (st.nservers > 1) os << "node{" << i << "}.setNumServers(" << fmt_d(st.nservers) << ");\n";
                break;
            }
            case NodeType::Router:
                os << "node{" << i << "} = Router(model, '" << nd.name << "');\n";
                break;
            case NodeType::Fork:
                os << "node{" << i << "} = Fork(model, '" << nd.name << "');\n";
                break;
            case NodeType::Join:
                os << "node{" << i << "} = Join(model, '" << nd.name << "', node{" << fork_of(sn, i) << "});\n";
                break;
            case NodeType::Sink:
                os << "node{" << i << "} = Sink(model, '" << nd.name << "');\n";
                break;
            case NodeType::ClassSwitch:
                os << "node{" << i << "} = Router(model, '" << nd.name
                   << "'); % Class switching is embedded in the routing matrix \n";
                break;
            default:
                break;
        }
    }
    os << "\n%% Block 2: classes\n";
    for (std::size_t k = 0; k < K; ++k) {
        const qn::JobClass& jc = sn.classes[k];
        const std::string kk = fmt_d(k + 1);
        if (std::isinf(jc.population)) {
            os << "jobclass{" << kk << "} = OpenClass(model, '" << jc.name << "', " << fmt_d(jc.prio) << ");\n";
        } else {
            const std::size_t ref = jc.population > 0 ? sn.station_to_node[jc.refstat - 1] : empty_class_ref(sn, k);
            os << "jobclass{" << kk << "} = ClosedClass(model, '" << jc.name << "', " << fmt_d(jc.population)
               << ", node{" << ref << "}, " << fmt_d(jc.prio) << ");\n";
        }
    }
    os << "\n";
    for (std::size_t i = 0; i < sn.stations.size(); ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const ProcSpec ps = proc_spec(sn, i, k);
            if (ps.kind == ProcSpec::SKIP) continue;
            const std::size_t nd = sn.station_to_node[i];
            const std::string head =
                "node{" + fmt_d(nd) + "}." + (ps.arrival ? "setArrival" : "setService") + "(jobclass{" + fmt_d(k + 1) + "}, ";
            const std::string tail = "); % (" + sn.nodes[nd - 1].name + "," + sn.classes[k].name + ")\n";
            switch (ps.kind) {
                case ProcSpec::REPLAYER: os << head << "Replayer('" << ps.file << "')" << tail; break;
                case ProcSpec::IMMEDIATE: os << head << "Immediate()" << tail; break;
                case ProcSpec::EXP: os << head << "Exp.fitMean(" << fmt_f(ps.mean) << ")" << tail; break;
                case ProcSpec::APH:
                    os << head << "APH.fitMeanAndSCV(" << fmt_f(ps.mean) << "," << fmt_f(ps.scv) << ")" << tail;
                    break;
                case ProcSpec::DISABLED: os << head << "Disabled.getInstance()" << tail; break;
                case ProcSpec::ERLANG:
                    os << head << "Erlang(" << fmt_f(ps.nphases / ps.mean) << "," << fmt_f(ps.nphases) << ")" << tail;
                    break;
                default: break;
            }
        }
    os << "\n%% Block 3: topology";
    os << "\n";
    os << "P = model.initRoutingMatrix(); % initialize routing matrix \n";
    for (const Route& r : routes(sn)) {
        const bool fork = sn.nodes[r.i].nodetype == NodeType::Fork;
        os << "P{" << r.k + 1 << "," << r.c + 1 << "}(" << r.i + 1 << "," << r.m + 1
           << ") = " << (fork ? std::string("1.0") : fmt_d(r.p)) << "; % (" << sn.nodes[r.i].name << ","
           << sn.classes[r.k].name << ") -> (" << sn.nodes[r.m].name << "," << sn.classes[r.c].name << ")\n";
    }
    os << "model.link(P);\n";
}

/** QN2MATLAB on a model under construction; refreshes it first. */
template <class T>
void qn2matlab(qn::Network<T>& model, const std::string& model_name, std::ostream& os) {
    qn2matlab(model.get_struct(), model_name, os);
}

/** QN2MATLAB with MATLAB's default model name. */
template <class T>
void qn2matlab(qn::Network<T>& model, std::ostream& os) {
    qn2matlab(model.get_struct(), "myModel", os);
}

/** QN2MATLAB into a file, as MATLAB does when `fid` is a file name. */
template <class T>
void qn2matlab(qn::Network<T>& model, const std::string& model_name, const std::string& path) {
    std::ofstream f = code_gen_detail::open_out(path);
    qn2matlab(model.get_struct(), model_name, f);
}

/** QN2MATLAB returned as a string. */
template <class T>
std::string qn2matlab_string(qn::Network<T>& model, const std::string& model_name = "myModel") {
    std::ostringstream os;
    qn2matlab(model.get_struct(), model_name, os);
    return os.str();
}

// ---------------------------------------------------------------------------
// QN2JAVA
// ---------------------------------------------------------------------------

/**
 * Port of `QN2JAVA(model, modelName, fid, headers)`: the body of a JLINE method
 * `public static Network ex()` that rebuilds the network. With `headers` false
 * only the statements are written, for pasting into an existing method.
 */
template <class T>
void qn2java(const qn::NetworkStruct<T>& sn, const std::string& model_name, std::ostream& os,
             bool headers = true) {
    using namespace code_gen_detail;
    using lang::NodeType;
    check_node_types(sn, "QN2JAVA");
    const std::size_t K = sn.classes.size();
    if (headers) os << "\tpublic static Network ex() {\n";
    os << "\t\tNetwork model = new Network(\"" << model_name << "\");\n";
    os << "\n\t\t// Block 1: nodes";
    os << "\t\t\t\n";
    for (std::size_t n = 0; n < sn.nodes.size(); ++n) {
        const std::size_t i = n + 1;
        const qn::NodeDef& nd = sn.nodes[n];
        switch (nd.nodetype) {
            case NodeType::Source:
                os << "\t\tSource node" << i << " = new Source(model, \"" << nd.name << "\");\n";
                break;
            case NodeType::Delay:
                os << "\t\tDelay node" << i << " = new Delay(model, \"" << nd.name << "\");\n";
                break;
            case NodeType::Queue: {
                const qn::Station<T>& st = sn.stations[nd.station - 1];
                os << "\t\tQueue node" << i << " = new Queue(model, \"" << nd.name << "\", SchedStrategy."
                   << sched_property(st.sched) << ");\n";
                if (st.nservers > 1) {
                    if (std::isinf(st.nservers))
                        os << "\t\tnode" << i << ".setNumberOfServers(Integer.MAX_VALUE);\n";
                    else
                        os << "\t\tnode" << i << ".setNumberOfServers(" << fmt_d(st.nservers) << ");\n";
                }
                break;
            }
            case NodeType::Router:
                os << "\t\tRouter node" << i << " = new Router(model, \"" << nd.name << "\");\n";
                break;
            case NodeType::Fork:
                os << "\t\tFork node" << i << " = new Fork(model, \"" << nd.name << "\");\n";
                break;
            case NodeType::Join:
                os << "\t\tJoin node" << i << " = new Join(model, \"" << nd.name << "\", node" << fork_of(sn, i)
                   << ");\n";
                break;
            case NodeType::Sink:
                os << "\t\tSink node" << i << " = new Sink(model, \"" << nd.name << "\");\n";
                break;
            case NodeType::ClassSwitch:
                os << "\t\tRouter node" << i << " = new Router(model, \"" << nd.name
                   << "\"); // Dummy node, class switching is embedded in the routing matrix P \n";
                break;
            default:
                break;
        }
    }
    os << "\n\t\t// Block 2: classes\n";
    for (std::size_t k = 0; k < K; ++k) {
        const qn::JobClass& jc = sn.classes[k];
        if (std::isinf(jc.population)) {
            os << "\t\tOpenClass jobclass" << k + 1 << " = new OpenClass(model, \"" << jc.name << "\", "
               << fmt_d(jc.prio) << ");\n";
        } else {
            const std::size_t ref = jc.population > 0 ? sn.station_to_node[jc.refstat - 1] : empty_class_ref(sn, k);
            os << "\t\tClosedClass jobclass" << k + 1 << " = new ClosedClass(model, \"" << jc.name << "\", "
               << fmt_d(jc.population) << ", node" << ref << ", " << fmt_d(jc.prio) << ");\n";
        }
    }
    os << "\t\t\n";
    for (std::size_t i = 0; i < sn.stations.size(); ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const ProcSpec ps = proc_spec(sn, i, k);
            if (ps.kind == ProcSpec::SKIP) continue;
            const std::size_t nd = sn.station_to_node[i];
            const std::string head = "\t\tnode" + fmt_d(nd) + "." + (ps.arrival ? "setArrival" : "setService") +
                                     "(jobclass" + fmt_d(k + 1) + ", ";
            const std::string tail = "); // (" + sn.nodes[nd - 1].name + "," + sn.classes[k].name + ")\n";
            // The weight argument is written only for a service, and not for an Immediate or Disabled one.
            double w = 1.0;
            const std::vector<T>& sp = sn.stations[i].schedparam;
            if (k < sp.size()) w = num_traits<T>::to_double(sp[k]);
            const std::string weight = (!ps.arrival && w != 1.0) ? ", " + fmt_f(w) : std::string();
            switch (ps.kind) {
                case ProcSpec::REPLAYER:
                    os << head << "new Replayer(\"" << ps.file << "\")" << weight << tail;
                    break;
                case ProcSpec::IMMEDIATE: os << head << "Immediate.getInstance()" << tail; break;
                case ProcSpec::EXP: os << head << "Exp.fitMean(" << fmt_f(ps.mean) << ")" << weight << tail; break;
                case ProcSpec::APH:
                    os << head << "APH.fitMeanAndSCV(" << fmt_f(ps.mean) << "," << fmt_f(ps.scv) << ")" << weight
                       << tail;
                    break;
                case ProcSpec::DISABLED: os << head << "Disabled.getInstance()" << tail; break;
                case ProcSpec::ERLANG:
                    os << head << "new Erlang(" << fmt_f(ps.nphases / ps.mean) << "," << fmt_d(ps.nphases) << ")"
                       << weight << tail;
                    break;
                default: break;
            }
        }
    os << "\n\t\t// Block 3: topology";
    os << "\t\n";
    os << "\t\tRoutingMatrix routingMatrix = model.initRoutingMatrix(); \n";
    os << "\t\n";
    for (const Route& r : routes(sn)) {
        const bool fork = sn.nodes[r.i].nodetype == NodeType::Fork;
        os << "\t\troutingMatrix.set(jobclass" << r.k + 1 << ", jobclass" << r.c + 1 << ", node" << r.i + 1
           << ", node" << r.m + 1 << ", " << fmt_f(fork ? sn.nodes[r.i].tasks_per_link : r.p) << "); // ("
           << sn.nodes[r.i].name << "," << sn.classes[r.k].name << ") -> (" << sn.nodes[r.m].name << ","
           << sn.classes[r.c].name << ")\n";
    }
    os << "\n\t\tmodel.link(routingMatrix);\n\n";
    if (headers) {
        os << "\t\treturn model;\n";
        os << "\t}\n";
    }
}

/** QN2JAVA on a model under construction; refreshes it first. */
template <class T>
void qn2java(qn::Network<T>& model, const std::string& model_name, std::ostream& os, bool headers = true) {
    qn2java(model.get_struct(), model_name, os, headers);
}

/** QN2JAVA with MATLAB's default model name. */
template <class T>
void qn2java(qn::Network<T>& model, std::ostream& os) {
    qn2java(model.get_struct(), "myModel", os, true);
}

/** QN2JAVA into a file. */
template <class T>
void qn2java(qn::Network<T>& model, const std::string& model_name, const std::string& path,
             bool headers = true) {
    std::ofstream f = code_gen_detail::open_out(path);
    qn2java(model.get_struct(), model_name, f, headers);
}

/** QN2JAVA returned as a string. */
template <class T>
std::string qn2java_string(qn::Network<T>& model, const std::string& model_name = "myModel",
                           bool headers = true) {
    std::ostringstream os;
    qn2java(model.get_struct(), model_name, os, headers);
    return os.str();
}

// ---------------------------------------------------------------------------
// LQN2JAVA
// ---------------------------------------------------------------------------

namespace code_gen_detail {

/** `jnum`: a Java double literal at the shortest spelling that reads back as `v`. */
inline std::string j_num(double v) {
    if (std::isnan(v)) return "Double.NaN";
    if (std::isinf(v)) return v > 0 ? "Double.POSITIVE_INFINITY" : "Double.NEGATIVE_INFINITY";
    char buf[40];
    if (v == std::round(v) && std::fabs(v) < 9007199254740992.0) {
        std::snprintf(buf, sizeof(buf), "%.0f.0", v == 0.0 ? 0.0 : v);
        return buf;
    }
    for (int p = 15; p <= 17; ++p) {
        std::snprintf(buf, sizeof(buf), "%.*g", p, v);
        if (std::strtod(buf, nullptr) == v) break;
    }
    return buf;
}

template <class T>
std::string j_num_t(const T& v) {
    return j_num(num_traits<T>::to_double(v));
}

/** `jarr`: a `new double[]{...}` literal. */
template <class T>
std::string j_arr(const std::vector<T>& v) {
    std::string s = "new double[]{";
    for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + j_num_t(v[i]);
    return s + "}";
}

/** `jmat` of a row vector: `new Matrix(new double[][]{{a, b}})`. */
template <class T>
std::string j_mat_row(const std::vector<T>& v) {
    std::string s = "new Matrix(new double[][]{{";
    for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + j_num_t(v[i]);
    return s + "}})";
}

/** `jmat` of a matrix, row by row. */
template <class T>
std::string j_mat(const Matrix<T>& m) {
    std::string s = "new Matrix(new double[][]{";
    for (std::size_t i = 0; i < m.rows(); ++i) {
        s += i ? ", {" : "{";
        for (std::size_t j = 0; j < m.cols(); ++j) s += (j ? ", " : "") + j_num_t(m(i, j));
        s += "}";
    }
    return s + "})";
}

inline std::string j_round(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.0f", std::round(v) == 0.0 ? 0.0 : std::round(v));
    return buf;
}

/**
 * `javaDist`: the Java constructor call rebuilding `d` with the same parameters.
 *
 * A law with no Java constructor listed here is written as an APH fitted to its
 * mean and SCV, where MATLAB also warns.
 */
template <class T>
std::string java_dist(const lang::Distrib<T>& d) {
    using lang::ProcessType;
    const std::vector<T>& p = d.params;
    auto par = [&](std::size_t k) -> const T& {
        if (k >= p.size())
            throw UnsupportedError(std::string("LQN2JAVA: the ") + lang::process_to_text(d.type) +
                                   " carries fewer parameters than its constructor takes");
        return p[k];
    };
    auto jn = [&](std::size_t k) { return j_num_t(par(k)); };
    auto jr = [&](std::size_t k) { return j_round(num_traits<T>::to_double(par(k))); };
    auto half = [&](std::size_t from, std::size_t n) { return std::vector<T>(p.begin() + from, p.begin() + from + n); };
    auto cls2 = [&](const char* c) { return std::string("new ") + c + "(" + jn(0) + ", " + jn(1) + ")"; };
    if (d.disabled || d.type == ProcessType::DISABLED) return "new Disabled()";
    switch (d.type) {
        case ProcessType::IMMEDIATE: return "new Immediate()";
        case ProcessType::EXP:
            return "new Exp(" + (p.empty() ? j_num(1.0 / num_traits<T>::to_double(d.mean)) : jn(0)) + ")";
        case ProcessType::DET: return "new Det(" + jn(0) + ")";
        case ProcessType::GEOMETRIC: return "new Geometric(" + jn(0) + ")";
        case ProcessType::POISSON: return "new Poisson(" + jn(0) + ")";
        case ProcessType::BERNOULLI: return "new Bernoulli(" + jn(0) + ")";
        case ProcessType::ERLANG: return "new Erlang(" + jn(0) + ", " + jr(1) + ")";
        case ProcessType::HYPEREXP:
            if (p.size() == 3) return "new HyperExp(" + jn(0) + ", " + jn(1) + ", " + jn(2) + ")";
            return "new HyperExp(" + j_arr(half(0, p.size() / 2)) + ", " + j_arr(half(p.size() / 2, p.size() / 2)) + ")";
        case ProcessType::GAMMA: return cls2("Gamma");
        case ProcessType::LOGNORMAL: return cls2("Lognormal");
        case ProcessType::UNIFORM: return cls2("Uniform");
        case ProcessType::PARETO: return cls2("Pareto");
        case ProcessType::NORMAL: return cls2("Normal");
        case ProcessType::DUNIFORM: return cls2("DiscreteUniform");
        case ProcessType::BINOMIAL: return "new Binomial(" + jr(0) + ", " + jn(1) + ")";
        case ProcessType::WEIBULL: return "new Weibull(" + jn(1) + ", " + jn(0) + ")";  // stored (scale, shape)
        case ProcessType::COXIAN:
        case ProcessType::COX2:
            // the 3-argument Coxian(mu1, mu2, phi1) is written as mu = [mu1, mu2], phi = [phi1, 1], which cox2() stores
            return "new Coxian(" + j_mat_row(half(0, p.size() / 2)) + ", " + j_mat_row(half(p.size() / 2, p.size() / 2)) + ")";
        case ProcessType::PH: return "new PH(" + j_mat_row(p) + ", " + j_mat(d.D0) + ")";
        case ProcessType::APH: return "new APH(" + j_mat_row(p) + ", " + j_mat(d.D0) + ")";
        case ProcessType::ME: return "new ME(" + j_mat_row(p) + ", " + j_mat(d.D0) + ")";
        case ProcessType::MAP: return "new MAP(" + j_mat(d.D0) + ", " + j_mat(d.D1) + ")";
        case ProcessType::RAP: return "new RAP(" + j_mat(d.D0) + ", " + j_mat(d.D1) + ")";
        case ProcessType::MMPP2: return "new MMPP2(" + jn(0) + ", " + jn(1) + ", " + jn(2) + ", " + jn(3) + ")";
        case ProcessType::ZIPF: return "new Zipf(" + jn(0) + ", " + jr(1) + ")";
        case ProcessType::DISCRETESAMPLER: {
            std::vector<T> x = d.trace;
            if (x.empty())
                for (std::size_t k = 0; k < p.size(); ++k) x.push_back(num_traits<T>::from_int(static_cast<long>(k + 1)));
            return "new DiscreteSampler(" + j_mat_row(p) + ", " + j_mat_row(x) + ")";
        }
        case ProcessType::REPLAYER: {
            if (d.trace_file.empty())
                throw UnsupportedError("LQN2JAVA: a Replayer built from in-memory samples has no trace file to name");
            std::string f;
            for (char c : d.trace_file) {
                if (c == '\\' || c == '"') f.push_back('\\');
                f.push_back(c);
            }
            return "new Replayer(\"" + f + "\")";
        }
        default:
            return "APH.fitMeanAndSCV(" + j_num_t(d.mean) + ", " + j_num_t(d.scv) + ")";
    }
}

/**
 * `sn.replygraph`: (nacts+1) x (nentries+1), 1-based, true where the activity
 * replies to the entry. Declared replies, plus the implicit leaf replies of
 * getStruct.m for an entry that declares none.
 */
template <class T>
std::vector<std::vector<bool>> reply_graph(const lqn::LqnModel<T>& m, const lqn::LqnStruct<T>& sn) {
    std::vector<std::vector<bool>> g(sn.nacts + 1, std::vector<bool>(sn.nentries + 1, false));
    std::map<std::string, std::size_t> ent, act;
    for (std::size_t e = 1; e <= sn.nentries; ++e) ent[sn.names[sn.eshift + e]] = e;
    for (std::size_t a = 1; a <= sn.nacts; ++a) act[sn.names[sn.ashift + a]] = a;
    const std::map<std::string, std::vector<std::string>> rep = lqn::detail::reply_activities(m, sn);
    for (const auto& kv : rep) {
        const auto ei = ent.find(kv.first);
        if (ei == ent.end()) continue;
        for (const std::string& an : kv.second) {
            const auto ai = act.find(an);
            if (ai != act.end()) g[ai->second][ei->second] = true;
        }
    }
    return g;
}

/** The quorum of the AND-join into `post_name` on task slot `tslot`, 0 when it declares none. */
template <class T>
std::size_t and_join_quorum(const lqn::LqnModel<T>& m, std::size_t tslot, const std::string& post_name) {
    if (tslot >= m.tasks.size()) return 0;
    for (const lqn::detail::RawPrecedence<T>& p : m.tasks[tslot].precedences) {
        if (p.pretype != lang::PrecedenceType::PRE_AND) continue;
        if (std::find(p.postacts.begin(), p.postacts.end(), post_name) == p.postacts.end()) continue;
        if (p.has_quorum) return p.quorum;
        if (!p.preparams.empty()) return static_cast<std::size_t>(num_traits<T>::to_double(p.preparams[0]));
        return 0;
    }
    return 0;
}

inline std::string upper_nospace(const std::string& s) {
    std::string out;
    for (std::size_t i = 0; i < s.size(); ++i)
        if (s[i] != ' ') out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(s[i]))));
    return out;
}

}  // namespace code_gen_detail

/**
 * Port of `LQN2JAVA(model, modelName, fid)`: a JLINE program that rebuilds the
 * layered network and solves it with SolverLN.
 *
 * It takes the intermediate `LqnModel` rather than the struct alone, because
 * the struct keeps neither the declared reply activities nor an AND-join's
 * quorum, both of which MATLAB reads off the model object.
 */
template <class T>
void lqn2java(const lqn::LqnModel<T>& model, const std::string& model_name, std::ostream& os) {
    using namespace code_gen_detail;
    using lang::PrecedenceType;
    const lqn::LqnStruct<T> sn = lqn::lqn_finalize(model);
    const std::vector<std::vector<bool>> replygraph = reply_graph(model, sn);
    const T zero = num_traits<T>::from_int(0);
    auto w = [&](std::size_t i, std::size_t j) { return num_traits<T>::to_double(sn.graph.get(i, j)); };

    os << "package jline.examples;\n\n";
    os << "import java.util.ArrayList;\n";
    os << "import jline.lang.*;\n";
    os << "import jline.lang.layered.*;\n";
    os << "import jline.lang.constant.*;\n";
    os << "import jline.lang.processes.*;\n";
    os << "import jline.util.matrix.Matrix;\n";
    os << "import jline.solvers.ln.SolverLN;\n\n";
    os << "public class TestSolver" << upper_nospace(model_name) << " {\n\n";
    os << "\tpublic static void main(String[] args) throws Exception{\n\n";
    os << "\tLayeredNetwork model = new LayeredNetwork(\"" << model_name << "\");\n";
    os << "\n";
    for (std::size_t h = 1; h <= sn.nhosts; ++h) {
        const std::string mult = std::isinf(sn.mult[h]) ? std::string("Integer.MAX_VALUE") : fmt_d(sn.mult[h]);
        os << "\tProcessor P" << h << " = new Processor(model, \"" << sn.names[h] << "\", " << mult << ", "
           << sched_feature(sn.sched[h]) << ");\n";
        if (sn.repl[h] != 1) os << "P" << h << ".setReplication(" << fmt_d(sn.repl[h]) << ");\n";
    }
    os << "\n";
    for (std::size_t t = 1; t <= sn.ntasks; ++t) {
        const std::size_t tidx = sn.tshift + t;
        const std::string mult =
            std::isinf(sn.mult[tidx]) ? std::string("Integer.MAX_VALUE") : fmt_d(sn.mult[tidx]);
        os << "\tTask T" << t << " = new Task(model, \"" << sn.names[tidx] << "\", " << mult << ", "
           << sched_feature(sn.sched[tidx]) << ").on(P" << sn.parent[tidx] << ");\n";
        if (sn.repl[tidx] != 1) os << "\tT" << t << ".setReplication(" << fmt_d(sn.repl[tidx]) << ");\n";
        const lang::Distrib<T>& th = sn.think[tidx];
        if (!th.disabled && th.type != lang::ProcessType::DISABLED)
            os << "\tT" << t << ".setThinkTime(" << java_dist(th) << ");\n";
    }
    os << "\n";
    for (std::size_t e = 1; e <= sn.nentries; ++e) {
        const std::size_t eidx = sn.eshift + e;
        os << "\tEntry E" << e << " = new Entry(model, \"" << sn.names[eidx] << "\").on(T"
           << sn.parent[eidx] - sn.tshift << ");\n";
    }
    os << "\n";
    for (std::size_t a = 1; a <= sn.nacts; ++a) {
        const std::size_t aidx = sn.ashift + a;
        const std::size_t tidx = sn.parent[aidx];
        std::string bound;
        for (std::size_t e = 1; e <= sn.nentries; ++e)
            if (sn.graph.get(sn.eshift + e, aidx) != zero) bound += ".boundTo(E" + std::to_string(e) + ")";
        std::string replies;
        if (sn.sched[tidx] != lang::SchedStrategy::REF) {
            std::vector<std::size_t> rt;
            bool all_nonref = true;
            for (std::size_t e = 1; e <= sn.nentries; ++e)
                if (replygraph[a][e]) {
                    rt.push_back(e);
                    if (sn.isref[sn.parent[sn.eshift + e]]) all_nonref = false;
                }
            if (!rt.empty() && all_nonref)
                for (std::size_t e : rt) replies += ".repliesTo(E" + std::to_string(e) + ")";
        }
        std::string calls;
        for (std::size_t c = 1; c <= sn.ncalls; ++c) {
            if (sn.callpair_src[c] != aidx) continue;
            const std::string target = std::to_string(sn.callpair_dst[c] - sn.eshift);
            const std::string mean = fmt_g(num_traits<T>::to_double(sn.callproc_mean[c]));
            if (sn.calltype[c] == lang::CallType::SYNC) calls += ".synchCall(E" + target + "," + mean + ")";
            else if (sn.calltype[c] == lang::CallType::ASYNC) calls += ".asynchCall(E" + target + "," + mean + ")";
        }
        os << "\tActivity A" << a << " = new Activity(model, \"" << sn.names[aidx] << "\", "
           << java_dist(sn.hostdem[aidx]) << ").on(T" << tidx - sn.tshift << ");";
        if (!bound.empty()) os << " A" << a << bound << ";";
        if (!calls.empty()) os << " A" << a << calls << ";";
        if (!replies.empty()) os << " A" << a << replies << ";";
        os << "\n";
        // activity think time; the Activity default is Immediate, so only a non-default one is written
        const lang::Distrib<T>& ath = sn.actthink[aidx];
        if (!ath.disabled && ath.type != lang::ProcessType::DISABLED && ath.type != lang::ProcessType::IMMEDIATE)
            os << "\tA" << a << ".setThinkTime(" << java_dist(ath) << ");\n";
    }
    os << "\n";

    // Sequential precedences
    for (std::size_t ai = 1; ai <= sn.nacts; ++ai) {
        const std::size_t aidx = sn.ashift + ai;
        const std::size_t tidx = sn.parent[aidx];
        for (std::size_t bidx : sn.graph.succ(aidx))
            if (bidx > sn.ashift && sn.actpretype[aidx] == PrecedenceType::PRE_SEQ &&
                sn.actposttype[bidx] == PrecedenceType::POST_SEQ)
                os << "\tT" << tidx - sn.tshift << ".addPrecedence(ActivityPrecedence.Serial(\"" << sn.names[aidx]
                   << "\", \"" << sn.names[bidx] << "\"));\n";
    }

    // Loop precedences (POST_LOOP): follow the chain of loop successors to the end activity, whose weight is 1/count
    bool has_pre_acts = false, has_post_acts = false;
    std::vector<bool> processed(sn.nacts + 1, false);
    for (std::size_t ai = 1; ai <= sn.nacts; ++ai) {
        const std::size_t aidx = sn.ashift + ai;
        const std::size_t tidx = sn.parent[aidx];
        if (processed[ai]) continue;
        for (std::size_t bidx : sn.graph.succ(aidx)) {
            if (!(bidx > sn.ashift && sn.actposttype[bidx] == PrecedenceType::POST_LOOP)) continue;
            if (processed[bidx - sn.ashift]) continue;
            const std::size_t loop_start = bidx;
            std::vector<std::string> names;
            std::size_t cur = loop_start;
            while (true) {
                names.push_back(sn.names[cur]);
                processed[cur - sn.ashift] = true;
                std::size_t end_idx = 0, next_idx = 0;
                for (std::size_t s : sn.graph.succ(cur)) {
                    if (s <= sn.ashift || sn.actposttype[s] != PrecedenceType::POST_LOOP) continue;
                    if (s == loop_start) continue;
                    const double wt = w(cur, s);
                    if (wt > 0 && wt < 1) end_idx = s;
                    else next_idx = s;
                }
                if (end_idx > 0) {
                    const double wt = w(cur, end_idx);
                    const double counts = wt > 0 ? 1.0 / wt : 1.0;
                    names.push_back(sn.names[end_idx]);
                    processed[end_idx - sn.ashift] = true;
                    os << "\n\t// Loop Activity Precedence \n";
                    if (!has_pre_acts) {
                        os << "\tArrayList<String> precActs = new ArrayList<String>();\n";
                        has_pre_acts = true;
                    } else {
                        os << "\tprecActs = new ArrayList<String>();\n";
                    }
                    for (const std::string& n : names) os << "\tprecActs.add(\"" << n << "\");\n";
                    os << "\tT" << tidx - sn.tshift << ".addPrecedence(ActivityPrecedence.Loop(\"" << sn.names[aidx]
                       << "\", precActs, Matrix.singleton(" << fmt_g(counts) << ")));\n";
                    break;
                } else if (next_idx > 0) {
                    cur = next_idx;
                } else {
                    break;
                }
            }
            break;
        }
    }

    // OrFork precedences (POST_OR)
    bool has_probs = false;
    {
        std::size_t prec_marker = 0;
        std::string prec_acts, prob_string;
        for (std::size_t ai = 1; ai <= sn.nacts; ++ai) {
            const std::size_t aidx = sn.ashift + ai;
            const std::size_t tidx = sn.parent[aidx];
            std::size_t prob_ctr = 0;
            for (std::size_t bidx : sn.graph.succ(aidx)) {
                if (!(bidx > sn.ashift && sn.actposttype[bidx] == PrecedenceType::POST_OR)) continue;
                const std::string add = "\tprecActs.add(\"" + sn.names[bidx] + "\");\n";
                ++prob_ctr;
                const std::string prob =
                    "\tprobs.set(0," + std::to_string(prob_ctr - 1) + "," + fmt_g(w(aidx, bidx)) + ");\n";
                if (prec_marker == 0) {
                    prec_marker = aidx - sn.ashift;
                    prec_acts = add;
                    prob_string = prob;
                } else {
                    prec_acts += add;
                    prob_string += prob;
                }
            }
            if (prec_marker > 0) {
                os << "\n\t// OrFork Activity Precedence \n";
                if (!has_pre_acts) {
                    os << "\tArrayList<String> precActs = new ArrayList<String>();\n";
                    has_pre_acts = true;
                } else {
                    os << "\tprecActs = new ArrayList<String>();\n";
                }
                if (!has_probs) {
                    os << "\tMatrix probs = new Matrix(1," << prob_ctr << ");\n";
                    has_probs = true;
                } else {
                    os << "\tprobs = new Matrix(1," << prob_ctr << ");\n";
                }
                os << prec_acts << prob_string;
                os << "\tT" << tidx - sn.tshift << ".addPrecedence(ActivityPrecedence.OrFork(\""
                   << sn.names[prec_marker + sn.ashift] << "\", precActs, probs));\n";
                prec_marker = 0;
            }
        }
    }

    // AndFork precedences (POST_AND)
    for (std::size_t ai = 1; ai <= sn.nacts; ++ai) {
        const std::size_t aidx = sn.ashift + ai;
        const std::size_t tidx = sn.parent[aidx];
        std::string post_acts;
        for (std::size_t bidx : sn.graph.succ(aidx))
            if (bidx > sn.ashift && sn.actposttype[bidx] == PrecedenceType::POST_AND)
                post_acts += "\tpostActs.add(\"" + sn.names[bidx] + "\");\n";
        if (post_acts.empty()) continue;
        os << "\n\t// AndFork Activity Precedence \n";
        if (!has_post_acts) {
            os << "\tArrayList<String> postActs = new ArrayList<String>();\n";
            has_post_acts = true;
        } else {
            os << "\t postActs = new ArrayList<String>();\n";
        }
        os << post_acts;
        os << "\tT" << tidx - sn.tshift << ".addPrecedence(ActivityPrecedence.AndFork(\"" << sn.names[aidx]
           << "\", postActs));\n";
    }

    // OrJoin (PRE_OR) and AndJoin (PRE_AND) precedences, scanned from the last activity back
    for (int pass = 0; pass < 2; ++pass) {
        const PrecedenceType want = pass == 0 ? PrecedenceType::PRE_OR : PrecedenceType::PRE_AND;
        for (std::size_t bi = sn.nacts; bi >= 1; --bi) {
            const std::size_t bidx = sn.ashift + bi;
            const std::size_t tidx = sn.parent[bidx];
            std::string prec_acts;
            for (std::size_t aidx : sn.graph.pred(bidx))
                if (aidx > sn.ashift && sn.actpretype[aidx] == want)
                    prec_acts += "\tprecActs.add(\"" + sn.names[aidx] + "\");\n";
            if (prec_acts.empty()) continue;
            os << (pass == 0 ? "\n\t// OrJoin Activity Precedence \n" : "\n\t// AndJoin Activity Precedence \n");
            if (!has_pre_acts) {
                os << "\tArrayList<String> precActs = new ArrayList<String>();\n";
                has_pre_acts = true;
            } else {
                os << "\tprecActs = new ArrayList<String>();\n";
            }
            os << prec_acts;
            const std::size_t local = tidx - sn.tshift;
            if (pass == 0) {
                os << "\tT" << local << ".addPrecedence(ActivityPrecedence.OrJoin(precActs, \"" << sn.names[bidx]
                   << "\"));\n";
            } else {
                const std::size_t q = and_join_quorum(model, local - 1, sn.names[bidx]);
                if (q == 0)
                    os << "\tT" << local << ".addPrecedence(ActivityPrecedence.AndJoin(precActs, \""
                       << sn.names[bidx] << "\"));\n";
                else
                    os << "\tT" << local << ".addPrecedence(ActivityPrecedence.AndJoin(precActs, \""
                       << sn.names[bidx] << "\", Matrix.singleton(" << j_round(double(q)) << ")));\n";
            }
        }
    }

    os << "\n\t// Model solution \n";
    os << "\tSolverLN solver = new SolverLN(model);\n";
    os << "\tsolver.getEnsembleAvg();\n";
    os << "\t}\n}\n";
}

/** LQN2JAVA on a model under construction. */
template <class T>
void lqn2java(const lqn::LqnBuilder<T>& b, const std::string& model_name, std::ostream& os) {
    lqn2java(b.model(), model_name, os);
}

/** LQN2JAVA with MATLAB's default model name. */
template <class T>
void lqn2java(const lqn::LqnModel<T>& model, std::ostream& os) {
    lqn2java(model, "myLayeredModel", os);
}

/** LQN2JAVA into a file. */
template <class T>
void lqn2java(const lqn::LqnModel<T>& model, const std::string& model_name, const std::string& path) {
    std::ofstream f = code_gen_detail::open_out(path);
    lqn2java(model, model_name, f);
}

/** LQN2JAVA returned as a string. */
template <class T>
std::string lqn2java_string(const lqn::LqnModel<T>& model, const std::string& model_name = "myLayeredModel") {
    std::ostringstream os;
    lqn2java(model, model_name, os);
    return os.str();
}

// ---------------------------------------------------------------------------
// LQN2MATLAB
// ---------------------------------------------------------------------------

namespace code_gen_detail {

/** `scalar2code`: a double at the shortest of 15..17 significant digits that reads back as itself. */
inline std::string m_scalar(double v) {
    if (std::isnan(v)) return "NaN";
    if (std::isinf(v)) return v > 0 ? "Inf" : "-Inf";
    if (v == std::round(v) && std::fabs(v) < 1e15) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.0f", v == 0.0 ? 0.0 : v);
        return buf;
    }
    char buf[40];
    for (int p = 15; p <= 17; ++p) {
        std::snprintf(buf, sizeof(buf), "%.*g", p, v);
        if (std::strtod(buf, nullptr) == v) break;
    }
    return buf;
}

template <class T>
std::string m_scalar_t(const T& v) {
    return m_scalar(num_traits<T>::to_double(v));
}

/** `num2code` of a row vector: `[a, b]`, a scalar when it has one element, `[]` when empty. */
template <class T>
std::string m_row(const std::vector<T>& v) {
    if (v.empty()) return "[]";
    if (v.size() == 1) return m_scalar_t(v[0]);
    std::string s = "[";
    for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + m_scalar_t(v[i]);
    return s + "]";
}

/** `num2code` of a matrix: `[a, b; c, d]`. */
template <class T>
std::string m_matrix(const Matrix<T>& M) {
    if (M.rows() == 0 || M.cols() == 0) return "[]";
    if (M.rows() == 1 && M.cols() == 1) return m_scalar_t(M(0, 0));
    std::string s = "[";
    for (std::size_t i = 0; i < M.rows(); ++i) {
        if (i) s += "; ";
        for (std::size_t j = 0; j < M.cols(); ++j) s += (j ? ", " : "") + m_scalar_t(M(i, j));
    }
    return s + "]";
}

/** `q(str)`: a MATLAB single-quoted literal. */
inline std::string m_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "''";
        else out.push_back(c);
    }
    return out + "'";
}

/** `cellstr2code`: `{'a', 'b'}`. */
inline std::string m_cellstr(const std::vector<std::string>& c) {
    std::string s = "{";
    for (std::size_t i = 0; i < c.size(); ++i) s += (i ? ", " : "") + m_quote(c[i]);
    return s + "}";
}

inline std::string m_repl(lang::ReplacementStrategy r) {
    switch (r) {
        case lang::ReplacementStrategy::RR: return "ReplacementStrategy.RR";
        case lang::ReplacementStrategy::FIFO: return "ReplacementStrategy.FIFO";
        case lang::ReplacementStrategy::SFIFO: return "ReplacementStrategy.SFIFO";
        case lang::ReplacementStrategy::LRU: return "ReplacementStrategy.LRU";
        case lang::ReplacementStrategy::HLRU: return "ReplacementStrategy.HLRU";
        case lang::ReplacementStrategy::CLIMB: return "ReplacementStrategy.CLIMB";
        case lang::ReplacementStrategy::QLRU: return "ReplacementStrategy.QLRU";
    }
    return "ReplacementStrategy.RR";
}

/**
 * `dist2code`: the constructor call that rebuilds `d` from its parameters, which
 * this port keeps in MATLAB getParam order. A law with no constructor spelling
 * here is written as the APH fit of its moments, as MATLAB does for a class it
 * does not know.
 */
template <class T>
std::string m_dist(const lang::Distrib<T>& d) {
    using lang::ProcessType;
    const std::vector<T>& p = d.params;
    auto par = [&](std::size_t k) -> std::string {
        if (k >= p.size())
            throw UnsupportedError(std::string("LQN2MATLAB: the ") + lang::process_to_text(d.type) +
                                   " carries fewer parameters than its constructor takes");
        return m_scalar_t(p[k]);
    };
    auto half = [&](std::size_t from, std::size_t n) {
        return m_row(std::vector<T>(p.begin() + from, p.begin() + from + n));
    };
    if (d.disabled || d.type == ProcessType::DISABLED) return "Disabled()";
    switch (d.type) {
        case ProcessType::IMMEDIATE: return "Immediate()";
        case ProcessType::EXP:
            return "Exp(" + (p.empty() ? m_scalar(1.0 / num_traits<T>::to_double(d.mean)) : par(0)) + ")";
        case ProcessType::DET: return "Det(" + par(0) + ")";
        case ProcessType::ERLANG: return "Erlang(" + par(0) + ", " + par(1) + ")";
        case ProcessType::HYPEREXP:
            if (p.size() == 3) return "HyperExp(" + par(0) + ", " + par(1) + ", " + par(2) + ")";
            return "HyperExp(" + half(0, p.size() / 2) + ", " + half(p.size() / 2, p.size() / 2) + ")";
        case ProcessType::GAMMA: return "Gamma(" + par(0) + ", " + par(1) + ")";
        case ProcessType::LOGNORMAL: return "Lognormal(" + par(0) + ", " + par(1) + ")";
        case ProcessType::UNIFORM: return "Uniform(" + par(0) + ", " + par(1) + ")";
        case ProcessType::PARETO: return "Pareto(" + par(0) + ", " + par(1) + ")";
        case ProcessType::NORMAL: return "Normal(" + par(0) + ", " + par(1) + ")";
        case ProcessType::BINOMIAL: return "Binomial(" + par(0) + ", " + par(1) + ")";
        case ProcessType::DUNIFORM: return "DiscreteUniform(" + par(0) + ", " + par(1) + ")";
        case ProcessType::WEIBULL: return "Weibull(" + par(1) + ", " + par(0) + ")";  // stored (scale, shape)
        case ProcessType::GEOMETRIC: return "Geometric(" + par(0) + ")";
        case ProcessType::POISSON: return "Poisson(" + par(0) + ")";
        case ProcessType::BERNOULLI: return "Bernoulli(" + par(0) + ")";
        case ProcessType::COXIAN:
        case ProcessType::COX2:
            // params are mu then phi; cox2() is MATLAB's 3-argument Coxian(mu1, mu2, phi1)
            if (d.cox_scalar_form && p.size() == 4) return "Coxian(" + par(0) + ", " + par(1) + ", " + par(2) + ")";
            return "Coxian(" + half(0, p.size() / 2) + ", " + half(p.size() / 2, p.size() / 2) + ")";
        case ProcessType::PH: return "PH(" + m_row(p) + ", " + m_matrix(d.D0) + ")";
        case ProcessType::APH: return "APH(" + m_row(p) + ", " + m_matrix(d.D0) + ")";
        case ProcessType::ME: return "ME(" + m_row(p) + ", " + m_matrix(d.D0) + ")";
        case ProcessType::MAP: return "MAP(" + m_matrix(d.D0) + ", " + m_matrix(d.D1) + ")";
        case ProcessType::RAP: return "RAP(" + m_matrix(d.D0) + ", " + m_matrix(d.D1) + ")";
        case ProcessType::MMPP2:
            return "MMPP2(" + par(0) + ", " + par(1) + ", " + par(2) + ", " + par(3) + ")";
        case ProcessType::ZIPF: return "Zipf(" + par(0) + ", " + par(1) + ")";
        case ProcessType::DISCRETESAMPLER: {
            std::vector<T> x = d.trace;
            if (x.empty())
                for (std::size_t k = 0; k < p.size(); ++k) x.push_back(num_traits<T>::from_int(static_cast<long>(k + 1)));
            return "DiscreteSampler(" + m_row(p) + ", " + m_row(x) + ")";
        }
        case ProcessType::REPLAYER:
            if (d.trace_file.empty())
                throw UnsupportedError("LQN2MATLAB: a Replayer built from in-memory samples has no trace file to name");
            return "Replayer(" + m_quote(d.trace_file) + ")";
        default:
            return "APH.fitMeanAndSCV(" + m_scalar_t(d.mean) + ", " + m_scalar_t(d.scv) + ")";
    }
}

/** True for a think / setup / delay-off time left at its constructor default. */
template <class T>
bool m_default_time(const lang::Distrib<T>& d) {
    return d.disabled || d.type == lang::ProcessType::DISABLED || d.type == lang::ProcessType::IMMEDIATE;
}

/** The `ActivityPrecedence` expression of one declared precedence (`precCode`). */
template <class T>
std::string m_prec(const lqn::detail::RawPrecedence<T>& ap) {
    using lang::PrecedenceType;
    const std::vector<std::string>& pre = ap.preacts;
    const std::vector<std::string>& post = ap.postacts;
    std::vector<T> preparams = ap.preparams;
    if (preparams.empty() && ap.has_quorum) preparams.push_back(num_traits<T>::from_int(static_cast<long>(ap.quorum)));
    const PrecedenceType a = ap.pretype, b = ap.posttype;
    if (a == PrecedenceType::PRE_SEQ && b == PrecedenceType::POST_SEQ && pre.size() == 1 && post.size() == 1 &&
        preparams.empty() && ap.postparams.empty())
        return "ActivityPrecedence.Serial(" + m_quote(pre[0]) + ", " + m_quote(post[0]) + ")";
    if (a == PrecedenceType::PRE_AND && b == PrecedenceType::POST_SEQ && post.size() == 1 && ap.postparams.empty()) {
        if (preparams.empty()) return "ActivityPrecedence.AndJoin(" + m_cellstr(pre) + ", " + m_quote(post[0]) + ")";
        return "ActivityPrecedence.AndJoin(" + m_cellstr(pre) + ", " + m_quote(post[0]) + ", " + m_row(preparams) + ")";
    }
    if (a == PrecedenceType::PRE_OR && b == PrecedenceType::POST_SEQ && post.size() == 1 && preparams.empty() &&
        ap.postparams.empty())
        return "ActivityPrecedence.OrJoin(" + m_cellstr(pre) + ", " + m_quote(post[0]) + ")";
    if (a == PrecedenceType::PRE_SEQ && b == PrecedenceType::POST_AND && pre.size() == 1 && preparams.empty() &&
        ap.postparams.empty())
        return "ActivityPrecedence.AndFork(" + m_quote(pre[0]) + ", " + m_cellstr(post) + ")";
    if (a == PrecedenceType::PRE_SEQ && b == PrecedenceType::POST_OR && pre.size() == 1 && preparams.empty())
        return "ActivityPrecedence.OrFork(" + m_quote(pre[0]) + ", " + m_cellstr(post) + ", " + m_row(ap.postparams) + ")";
    if (a == PrecedenceType::PRE_SEQ && b == PrecedenceType::POST_LOOP && pre.size() == 1 && preparams.empty()) {
        // The loop count is ONE number; this port stores it once per body activity.
        std::vector<T> count;
        if (!ap.postparams.empty()) count.push_back(ap.postparams[0]);
        return "ActivityPrecedence.Loop(" + m_quote(pre[0]) + ", " + m_cellstr(post) + ", " + m_row(count) + ")";
    }
    if (a == PrecedenceType::PRE_SEQ && b == PrecedenceType::POST_CACHE && pre.size() == 1 && preparams.empty() &&
        ap.postparams.empty())
        return "ActivityPrecedence.CacheAccess(" + m_quote(pre[0]) + ", " + m_cellstr(post) + ")";
    auto type_name = [](PrecedenceType t) -> std::string {
        switch (t) {
            case PrecedenceType::PRE_SEQ: return "ActivityPrecedenceType.PRE_SEQ";
            case PrecedenceType::PRE_AND: return "ActivityPrecedenceType.PRE_AND";
            case PrecedenceType::PRE_OR: return "ActivityPrecedenceType.PRE_OR";
            case PrecedenceType::POST_SEQ: return "ActivityPrecedenceType.POST_SEQ";
            case PrecedenceType::POST_AND: return "ActivityPrecedenceType.POST_AND";
            case PrecedenceType::POST_OR: return "ActivityPrecedenceType.POST_OR";
            case PrecedenceType::POST_LOOP: return "ActivityPrecedenceType.POST_LOOP";
            case PrecedenceType::POST_CACHE: return "ActivityPrecedenceType.POST_CACHE";
        }
        return m_scalar(double(static_cast<int>(t)));
    };
    return "ActivityPrecedence(" + m_cellstr(pre) + ", " + m_cellstr(post) + ", " + type_name(a) + ", " +
           type_name(b) + ", " + m_row(preparams) + ", " + m_row(ap.postparams) + ")";
}

/** `emitServerExtras`: admission constraints, load dependence and server pools of a host or task. */
template <class T>
void m_server_extras(std::ostream& os, const std::string& v, const std::string& owner, const Matrix<T>* A,
                     const std::vector<T>* b, const std::vector<lqn::detail::RawLinConRow<T>>* rows,
                     const std::vector<T>* lld, bool has_cd, bool has_jd,
                     const std::vector<lqn::detail::RawServerPool<T>>* pools) {
    if (A && A->rows() > 0 && A->cols() > 0) os << v << ".setConstraint(" << m_matrix(*A) << ", " << m_row(*b) << ");\n";
    if (rows)
        for (const lqn::detail::RawLinConRow<T>& r : *rows)
            os << v << ".addConstraint(" << m_cellstr(r.names) << ", " << m_row(r.coeffs) << ", " << m_scalar_t(r.cap)
               << ");\n";
    if (lld && !lld->empty()) os << v << ".setLoadDependence(" << m_row(*lld) << ");\n";
    if (has_cd || has_jd)
        throw UnsupportedError("LQN2MATLAB: " + owner + " declares a " + (has_cd ? "class" : "joint") +
                               "-dependent rate as a compiled function, which has no MATLAB source to write");
    if (pools)
        for (const lqn::detail::RawServerPool<T>& sp : *pools)
            os << v << ".addServerType(ServerType(" << m_quote(sp.name) << ", " << m_scalar(sp.count) << ", "
               << m_cellstr(sp.compatible) << ", " << m_scalar_t(sp.rate) << "));\n";
}

}  // namespace code_gen_detail

/**
 * Port of `LQN2MATLAB(model, modelName, fid)`, the generator form: a MATLAB
 * script that rebuilds the layered network by the constructors a user writes,
 * with hosts, tasks, entries and activities in model order so the regenerated
 * model has the same LayeredNetworkStruct indexing, and numbers at the shortest
 * spelling that reads back bit for bit.
 *
 * It reads the intermediate `LqnModel`, which is the C++ counterpart of the
 * MATLAB objects. What that model does not carry is written at its MATLAB
 * default: a task priority, an entry type other than PH1PH2, an activity call
 * order, phase, a processor class other than Processor, and the numeric-versus-
 * object choice MATLAB makes for an Immediate or Exp time (this port writes the
 * object form, and cannot tell an explicit `Immediate()` think time from the
 * default one). A class- or joint-dependent rate is a compiled function here
 * and is refused.
 */
template <class T>
void lqn2matlab(const lqn::LqnModel<T>& m, const std::string& model_name, std::ostream& os) {
    using namespace code_gen_detail;
    std::map<std::string, std::size_t> entry_idx, act_idx;
    for (std::size_t e = 0; e < m.entries.size(); ++e) entry_idx.emplace(m.entries[e].name, e + 1);
    for (std::size_t a = 0; a < m.acts.size(); ++a) act_idx.emplace(m.acts[a].name, a + 1);
    auto entry_ref = [&](const std::string& n) {
        const auto it = entry_idx.find(n);
        return it == entry_idx.end() ? m_quote(n) : "E{" + std::to_string(it->second) + "}";
    };
    auto sched = [](lang::SchedStrategy s) { return "SchedStrategy." + sched_property(s); };

    os << "% LayeredNetwork generated by LQN2MATLAB\n";
    os << "model = LayeredNetwork(" << m_quote(model_name) << ");\n";
    os << "P = {}; T = {}; E = {}; A = {};\n";

    os << "\n%% Block 1: processors\n";
    for (std::size_t p = 0; p < m.procs.size(); ++p) {
        const lqn::detail::RawProc& h = m.procs[p];
        const std::string hv = "P{" + std::to_string(p + 1) + "}";
        // A quantum of 0 is this port's "not declared" (the .lqnx reader's default); MATLAB's default is 0.001.
        const double quantum = h.quantum == 0.0 ? 0.001 : h.quantum;
        const char* ctor = h.is_host_class ? " = Host(model, " : " = Processor(model, ";
        if (quantum != 0.001 || h.speed_factor != 1.0)
            os << hv << ctor << m_quote(h.name) << ", " << m_scalar(h.mult) << ", " << sched(h.sched)
               << ", " << m_scalar(quantum) << ", " << m_scalar(h.speed_factor) << ");\n";
        else
            os << hv << ctor << m_quote(h.name) << ", " << m_scalar(h.mult) << ", " << sched(h.sched)
               << ");\n";
        if (h.repl != 1.0) os << hv << ".setReplication(" << m_scalar(h.repl) << ");\n";
        auto lc = m.proc_lincon.find(p);
        auto rows = m.proc_linconrows.find(p);
        auto lld = m.proc_lldscaling.find(p);
        auto pools = m.proc_pools.find(p);
        auto cd = m.proc_cdscaling.find(p);
        auto jd = m.proc_jdscaling.find(p);
        m_server_extras<T>(os, hv, h.name, lc == m.proc_lincon.end() ? nullptr : &lc->second.first,
                           lc == m.proc_lincon.end() ? nullptr : &lc->second.second,
                           rows == m.proc_linconrows.end() ? nullptr : &rows->second,
                           lld == m.proc_lldscaling.end() ? nullptr : &lld->second,
                           cd != m.proc_cdscaling.end() && static_cast<bool>(cd->second),
                           jd != m.proc_jdscaling.end() && static_cast<bool>(jd->second),
                           pools == m.proc_pools.end() ? nullptr : &pools->second);
    }

    os << "\n%% Block 2: tasks\n";
    for (std::size_t t = 0; t < m.tasks.size(); ++t) {
        const lqn::detail::RawTask<T>& tk = m.tasks[t];
        const std::string tv = "T{" + std::to_string(t + 1) + "}";
        const std::string on = ".on(P{" + std::to_string(tk.proc_slot + 1) + "});\n";
        const bool setup = !m_default_time(tk.setuptime) || !m_default_time(tk.delayofftime);
        if (tk.nitems > 0) {
            std::vector<double> cap(tk.itemcap.begin(), tk.itemcap.end());
            os << tv << " = CacheTask(model, " << m_quote(tk.name) << ", " << m_scalar(double(tk.nitems)) << ", "
               << m_row(cap) << ", " << m_repl(tk.replacestrat) << ", " << m_scalar(tk.mult) << ", " << sched(tk.sched)
               << ")" << on;
            if (tk.retrieval) os << tv << ".setRetrieval(true);\n";
        } else {
            os << tv << " = " << (setup ? "SetupTask" : "Task") << "(model, " << m_quote(tk.name) << ", "
               << m_scalar(tk.mult) << ", " << sched(tk.sched) << ")" << on;
        }
        if (!m_default_time(tk.thinktime)) os << tv << ".setThinkTime(" << m_dist(tk.thinktime) << ");\n";
        if (!m_default_time(tk.setuptime)) os << tv << ".setSetupTime(" << m_dist(tk.setuptime) << ");\n";
        if (!m_default_time(tk.delayofftime)) os << tv << ".setDelayOffTime(" << m_dist(tk.delayofftime) << ");\n";
        if (tk.repl != 1.0) os << tv << ".setReplication(" << m_scalar(tk.repl) << ");\n";
        if (tk.priority != 0) os << tv << ".setPriority(" << tk.priority << ");\n";
        for (const auto& fi : tk.fanin) os << tv << ".setFanIn(" << m_quote(fi.first) << ", " << m_scalar(fi.second) << ");\n";
        for (const auto& fo : tk.fanout)
            os << tv << ".setFanOut(" << m_quote(fo.first) << ", " << m_scalar(fo.second) << ");\n";
        m_server_extras<T>(os, tv, tk.name, &tk.lincon_A, &tk.lincon_b, &tk.linconrows, &tk.lldscaling,
                           static_cast<bool>(tk.cdscaling), static_cast<bool>(tk.jdscaling), &tk.pools);
    }

    os << "\n%% Block 3: entries\n";
    for (std::size_t e = 0; e < m.entries.size(); ++e) {
        const lqn::detail::RawEntry<T>& en = m.entries[e];
        const std::string ev = "E{" + std::to_string(e + 1) + "}";
        const std::string on = ".on(T{" + std::to_string(en.task_slot + 1) + "});\n";
        if (en.cardinality > 0) {
            std::vector<T> x;
            for (std::size_t k = 0; k < en.popularity.size(); ++k) x.push_back(num_traits<T>::from_int(static_cast<long>(k + 1)));
            os << ev << " = ItemEntry(model, " << m_quote(en.name) << ", " << m_scalar(double(en.cardinality))
               << ", DiscreteSampler(" << m_row(en.popularity) << ", " << m_row(x) << "))" << on;
        } else {
            os << ev << " = Entry(model, " << m_quote(en.name) << ")" << on;
        }
        if (en.type != "PH1PH2") os << ev << ".setType(" << m_quote(en.type) << ");\n";
        if (en.has_arrival) os << ev << ".setArrival(" << m_dist(en.arrival) << ");\n";
    }
    for (std::size_t e = 0; e < m.entries.size(); ++e)
        for (std::size_t f = 0; f < m.entries[e].fwd_dest.size(); ++f)
            os << "E{" << e + 1 << "}.forward(" << entry_ref(m.entries[e].fwd_dest[f]) << ", "
               << m_scalar_t(m.entries[e].fwd_prob[f]) << ");\n";

    os << "\n%% Block 4: activities\n";
    for (std::size_t a = 0; a < m.acts.size(); ++a) {
        const lqn::detail::RawActivity<T>& ac = m.acts[a];
        const std::string av = "A{" + std::to_string(a + 1) + "}";
        os << av << " = Activity(model, " << m_quote(ac.name) << ", " << m_dist(ac.hostdem) << ").on(T{"
           << ac.task_slot + 1 << "})";
        if (!ac.bound_to_entry.empty()) os << ".boundTo(" << entry_ref(ac.bound_to_entry) << ")";
        os << ";\n";
        if (ac.call_order != "STOCHASTIC") os << av << ".setCallOrder(" << m_quote(ac.call_order) << ");\n";
        if (!m_default_time(ac.thinktime)) os << av << ".setThinkTime(" << m_dist(ac.thinktime) << ");\n";
        if (ac.phase != 1) os << av << ".setPhase(" << ac.phase << ");\n";
        for (const lqn::detail::RawCall<T>& c : ac.sync_calls)
            os << av << ".synchCall(" << entry_ref(c.dest) << ", " << m_scalar_t(c.mean) << ");\n";
        for (const auto& g : ac.call_groups) {
            std::string rs;
            if (g.first == lang::RoutingStrategy::RROBIN) rs = "RoutingStrategy.RROBIN";
            else if (g.first == lang::RoutingStrategy::JSQ) rs = "RoutingStrategy.JSQ";
            else
                throw UnsupportedError(std::string("LQN2MATLAB: call groups carry RROBIN or JSQ; routing strategy ") +
                                       lang::routing_to_text(g.first) + " cannot be written out");
            os << av << ".recordCallGroup(" << rs << ", " << m_cellstr(g.second) << ");\n";
        }
        for (const lqn::detail::RawCall<T>& c : ac.async_calls)
            os << av << ".asynchCall(" << entry_ref(c.dest) << ", " << m_scalar_t(c.mean) << ");\n";
    }

    os << "\n%% Block 5: replies\n";
    for (std::size_t e = 0; e < m.entries.size(); ++e)
        for (const std::string& rn : m.entries[e].reply_activities) {
            const auto it = act_idx.find(rn);
            const bool is_ref = it != act_idx.end() &&
                                m.tasks[m.acts[it->second - 1].task_slot].sched == lang::SchedStrategy::REF;
            if (it == act_idx.end() || is_ref)
                os << "E{" << e + 1 << "}.replyActivity{end+1} = " << m_quote(rn) << ";\n";
            else
                os << "A{" << it->second << "}.repliesTo(E{" << e + 1 << "});\n";
        }

    os << "\n%% Block 6: precedences\n";
    for (std::size_t t = 0; t < m.tasks.size(); ++t)
        for (const lqn::detail::RawPrecedence<T>& ap : m.tasks[t].precedences)
            os << "T{" << t + 1 << "}.addPrecedence(" << m_prec(ap) << ");\n";
}

/** LQN2MATLAB under the model's own name (`myLayeredModel` when it has none). */
template <class T>
void lqn2matlab(const lqn::LqnModel<T>& model, std::ostream& os) {
    lqn2matlab(model, model.name.empty() ? std::string("myLayeredModel") : model.name, os);
}

/** LQN2MATLAB on a model under construction. */
template <class T>
void lqn2matlab(const lqn::LqnBuilder<T>& b, const std::string& model_name, std::ostream& os) {
    lqn2matlab(b.model(), model_name, os);
}

/** LQN2MATLAB into a file. */
template <class T>
void lqn2matlab(const lqn::LqnModel<T>& model, const std::string& model_name, const std::string& path) {
    std::ofstream f = code_gen_detail::open_out(path);
    lqn2matlab(model, model_name, f);
}

/** LQN2MATLAB returned as a string. */
template <class T>
std::string lqn2matlab_string(const lqn::LqnModel<T>& model, const std::string& model_name) {
    std::ostringstream os;
    lqn2matlab(model, model_name, os);
    return os.str();
}

// ---------------------------------------------------------------------------
// LINE2MATLAB / LINE2JAVA
// ---------------------------------------------------------------------------

/** Port of `LINE2MATLAB(model)` for a Network: QN2MATLAB under the model's own name. */
template <class T>
void line2matlab(qn::Network<T>& model, std::ostream& os) {
    const qn::NetworkStruct<T>& sn = model.get_struct();
    qn2matlab(sn, sn.name, os);
}

/** Port of `LINE2MATLAB(model, filename)` for a Network. */
template <class T>
void line2matlab(qn::Network<T>& model, const std::string& path) {
    std::ofstream f = code_gen_detail::open_out(path);
    line2matlab(model, f);
}

/** Port of `LINE2JAVA(model)` for a Network: QN2JAVA under the model's own name. */
template <class T>
void line2java(qn::Network<T>& model, std::ostream& os) {
    const qn::NetworkStruct<T>& sn = model.get_struct();
    qn2java(sn, sn.name, os, true);
}

/** Port of `LINE2JAVA(model, filename)` for a Network. */
template <class T>
void line2java(qn::Network<T>& model, const std::string& path) {
    std::ofstream f = code_gen_detail::open_out(path);
    line2java(model, f);
}

/**
 * Port of `LINE2JAVA(model)` for a LayeredNetwork. The intermediate model does
 * not carry the model name, so it is passed explicitly.
 */
template <class T>
void line2java(const lqn::LqnModel<T>& model, const std::string& model_name, std::ostream& os) {
    lqn2java(model, model_name, os);
}

/** Port of `LINE2JAVA(model, filename)` for a LayeredNetwork. */
template <class T>
void line2java(const lqn::LqnModel<T>& model, const std::string& model_name, const std::string& path) {
    lqn2java(model, model_name, path);
}

/** Port of `LINE2MATLAB(model)` for a LayeredNetwork: LQN2MATLAB under the model's own name. */
template <class T>
void line2matlab(const lqn::LqnModel<T>& model, std::ostream& os) {
    lqn2matlab(model, os);
}

/** Port of `LINE2MATLAB(model, filename)` for a LayeredNetwork. */
template <class T>
void line2matlab(const lqn::LqnModel<T>& model, const std::string& path) {
    std::ofstream f = code_gen_detail::open_out(path);
    lqn2matlab(model, f);
}

/** Port of `LINE2JAVA(model)` for a LayeredNetwork under its own name (`myLayeredModel` when it has none). */
template <class T>
void line2java(const lqn::LqnModel<T>& model, std::ostream& os) {
    lqn2java(model, model.name.empty() ? std::string("myLayeredModel") : model.name, os);
}

namespace code_gen_detail {

/** The `name` of the model in a model.json envelope, `fallback` when it has none. */
inline std::string json_model_name(const std::string& path, const std::string& fallback) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("code_gen: cannot open " + path);
    detail::json root;
    try {
        in >> root;
    } catch (const detail::json::parse_error& e) {
        throw InputError("code_gen: malformed JSON in " + path + ": " + e.what());
    }
    const detail::json& model = root.contains("model") ? root.at("model") : root;
    return model.contains("name") && model.at("name").is_string() ? model.at("name").get<std::string>()
                                                                   : fallback;
}

}  // namespace code_gen_detail

/**
 * LINE2JAVA on a model.json file, dispatching on the model type as MATLAB
 * dispatches on the class of the object: a LayeredNetwork goes to LQN2JAVA and
 * any other network to QN2JAVA, each under the name the document declares.
 */
template <class T = double>
void line2java_json(const std::string& json_path, std::ostream& os) {
    if (is_layered_json(json_path)) {
        lqn2java(read_lqn_json_model<T>(json_path), code_gen_detail::json_model_name(json_path, "myLayeredModel"),
                 os);
    } else {
        qn::Network<T> model = read_network_json<T>(json_path);
        line2java(model, os);
    }
}

/** LINE2MATLAB on a model.json file: a LayeredNetwork goes to LQN2MATLAB, any other network to QN2MATLAB. */
template <class T = double>
void line2matlab_json(const std::string& json_path, std::ostream& os) {
    if (is_layered_json(json_path)) {
        lqn2matlab(read_lqn_json_model<T>(json_path), os);
        return;
    }
    qn::Network<T> model = read_network_json<T>(json_path);
    line2matlab(model, os);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_CODE_GEN_H
