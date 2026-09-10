/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_JMT_WRITER_H
#define LINE_IO_JMT_WRITER_H

/**
 * Port of `@@JMTIO`: a refreshed `NetworkStruct` written out as a JMT `.jsimg`
 * simulation model.
 *
 * THE REFERENCE WALKS NODE OBJECTS; THIS PORT WALKS THE STRUCT. MATLAB's
 * `writeJSIM` iterates `model.nodes{i}` and asks each for its three SECTIONS
 * (input, server, output), which are objects the node's constructor installed.
 * This port has no node objects -- `network_builder` produces a `NetworkStruct`
 * directly -- so `jmt_sections` re-derives the same triple from the node type
 * and the scheduling strategy, transcribing the constructors of `Queue.m`,
 * `Source.m`, `Sink.m`, `Router.m`, `ClassSwitch.m`, `Cache.m`, `Logger.m`,
 * `Fork.m`, `Join.m`, `Place.m` and `Transition.m`. That table is the one place
 * where this port can drift from the reference without a compiler error, so it
 * names its source file for each row.
 *
 * WHAT IS DELIBERATELY REPRODUCED RATHER THAN CORRECTED. Two places in the
 * reference emit less than they appear to, and all four codebases agree on the
 * result, so this port agrees with them and says so at the site:
 *   - `saveForkStrategy` emits ONE OutPath entry (the last connected node),
 *     not one per branch -- harmless only because `isSimplifiedFork` is true.
 *   - the analytic `MMPP2Par` branch of `saveServiceStrategy` is unreachable,
 *     because the MAP branch catches MMPP2 first.
 * Changing either here would make the C++ row the odd one out, which is the
 * opposite of what a port is for.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "line/io/jmt_dist.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/xml.h"

namespace line {
namespace io {

using lang::JoinStrategy;
using lang::NodeType;
using lang::PollingType;
using lang::RoutingStrategy;
using lang::SchedStrategy;
using qn::DropStrategy;

/**
 * The simulation controls the JSIM header carries, MATLAB's `JMTIO` properties.
 *
 * `max_events = -1` and `sim_conf_int = 0.99` / `sim_max_rel_err = 0.03` are the
 * reference's constructor defaults; `max_samples` is `options.samples`, which
 * `runAnalyzer` raises to 5000 before it reaches here.
 */
struct JmtWriteOptions {
    std::string file_name = "model";  ///< base name; the header echoes it plus `.jsimg`
    std::string log_path;             ///< `model.getLogPath`, the `logPath` attribute
    long seed = 23000;
    double max_samples = 10000.0;
    double max_events = -1.0;
    double max_simulated_time = std::numeric_limits<double>::infinity();
    double sim_conf_int = 0.99;
    double sim_max_rel_err = 0.03;
};

/**
 * `sn.connmatrix`: which ordered pairs of nodes the model LINKED.
 *
 * The reference carries the matrix on the struct, recorded by `Network.addLink`.
 * This port records the links only through the routing blocks `P`, so the
 * connection set is their union support. That is the same set for every model
 * built the way `link(P)` builds one -- `addLink` is called from `link` for each
 * nonzero entry -- and it is what the JSIM `<connection>` list, the ClassSwitch
 * row sum, the WRROBIN and PROB destination lists and the SPN input/output
 * vectors all read.
 *
 * `Peff` IS NOT USED HERE, deliberately: it has the RAND/RROBIN expansion and
 * the class-switch folding applied, so a Router's declared links would come
 * back as the links of the stations it routes to, and the exported topology
 * would not be the user's.
 */
template <class T>
std::vector<std::vector<bool>> jmt_conn_matrix(const qn::NetworkStruct<T>& sn) {
    const std::size_t I = sn.nodes.size();
    std::vector<std::vector<bool>> C(I, std::vector<bool>(I, false));
    const T zero = num_traits<T>::from_int(0);
    for (const auto& kv : sn.P) {
        const Matrix<T>& B = kv.second;
        for (std::size_t a = 0; a < B.rows() && a < I; ++a)
            for (std::size_t b = 0; b < B.cols() && b < I; ++b)
                if (!(B(a, b) == zero)) C[a][b] = true;
    }
    return C;
}

/**
 * Port of `getExportableClasses`.
 *
 * A closed class with no customers is dropped from the exported model UNLESS
 * jobs can still reach it: as a cache hit/miss class, or through a class switch
 * from a populated class of the same chain. Exporting it anyway would give JMT
 * a class it can never see and, for a closed class, a reference station with
 * zero population -- which JMT reports as an unreachable measure rather than as
 * an empty one.
 */
template <class T>
std::vector<bool> jmt_exportable_classes(const qn::NetworkStruct<T>& sn) {
    const std::size_t K = sn.nclasses;
    std::set<std::size_t> cache_classes;  // 1-based
    for (const auto& kv : sn.nodeparam) {
        const qn::CacheParam<T>& cp = kv.second;
        for (std::size_t r = 0; r < cp.hitclass.size(); ++r)
            if (cp.hitclass[r] > 0) cache_classes.insert(cp.hitclass[r]);
        for (std::size_t r = 0; r < cp.missclass.size(); ++r)
            if (cp.missclass[r] > 0) cache_classes.insert(cp.missclass[r]);
    }
    std::set<std::size_t> switch_classes;
    for (std::size_t c = 0; c < sn.chains.size(); ++c) {
        std::vector<std::size_t> in;
        for (std::size_t r = 0; r < K; ++r)
            if (sn.chains[c][r]) in.push_back(r + 1);
        if (in.size() <= 1) continue;
        bool has_jobs = false;
        for (std::size_t r : in)
            if (sn.classes[r - 1].population > 0.0) has_jobs = true;
        if (has_jobs)
            for (std::size_t r : in) switch_classes.insert(r);
    }
    std::vector<bool> keep(K, true);
    for (std::size_t r = 1; r <= K; ++r) {
        const double n = sn.classes[r - 1].population;
        if (std::isfinite(n) && n == 0.0 && !cache_classes.count(r) && !switch_classes.count(r))
            keep[r - 1] = false;
    }
    return keep;
}

/**
 * The three JMT section class names of a node: input, server, output.
 *
 * An empty entry is a section the node does not have -- only a Sink, whose
 * `getSections` returns `{'', JobSink, ''}`. The names are LINE's section class
 * names; `jmt_write_jsim` overwrites the ones JMT spells differently, exactly
 * where the reference does.
 */
template <class T>
struct JmtSections {
    std::string input, server, output;
};

template <class T>
JmtSections<T> jmt_sections(const qn::NetworkStruct<T>& sn, std::size_t ind) {
    JmtSections<T> s;
    const qn::NodeDef& nd = sn.nodes[ind - 1];
    const std::size_t ist = nd.station;
    switch (nd.nodetype) {
        case NodeType::Source:  // Source.m:66-68
            s.input = "RandomSource";
            s.server = "ServiceTunnel";
            s.output = "Dispatcher";
            break;
        case NodeType::Sink:  // Sink.m:28, getSections
            s.server = "JobSink";
            break;
        case NodeType::Router:  // Router.m:62-66
            s.input = "Buffer";
            s.server = "ServiceTunnel";
            s.output = "Dispatcher";
            break;
        case NodeType::ClassSwitch:  // ClassSwitch.m:32-38
            s.input = "Buffer";
            s.server = "StatelessClassSwitcher";
            s.output = "Dispatcher";
            break;
        case NodeType::Cache:  // Cache.m:57-58; the server section is a
            s.input = "Buffer";  // CacheClassSwitcher whose className is 'Cache'
            s.server = "Cache";
            s.output = "Dispatcher";
            break;
        case NodeType::Logger:  // Logger.m:40-47
            s.input = "Buffer";
            s.server = "LogTunnel";
            s.output = "Dispatcher";
            break;
        case NodeType::Fork:  // Fork.m:106-109
            s.input = "Buffer";
            s.server = "ServiceTunnel";
            s.output = "Forker";
            break;
        case NodeType::Join:  // Join.m:34-36
            s.input = "Joiner";
            s.server = "ServiceTunnel";
            s.output = "Dispatcher";
            break;
        case NodeType::Transition:  // Transition.m:29-36
            s.input = "Enabling";
            s.server = "Timing";
            s.output = "Firing";
            break;
        case NodeType::Place:  // Place.m:27-31 and installQueueServer
            s.input = "Storage";
            s.output = "Linkage";
            s.server = "ServiceTunnel";
            if (ist != 0) {
                bool queueing = false;
                for (std::size_t r = 0; r < sn.nclasses; ++r)
                    if (!sn.service[ist - 1][r].disabled) queueing = true;
                if (queueing) {
                    const SchedStrategy sch = sn.stations[ist - 1].sched;
                    if (sch == SchedStrategy::INF)
                        s.server = "InfiniteServer";
                    else if (sch == SchedStrategy::PS || sch == SchedStrategy::DPS ||
                             sch == SchedStrategy::GPS || sch == SchedStrategy::LPS)
                        s.server = "SharedServer";
                    else
                        s.server = "Server";
                }
            }
            break;
        case NodeType::Queue:
        case NodeType::Delay:
        default: {  // Queue.m:48-53 and the schedStrategy switch at :79-110
            s.input = "Buffer";
            s.output = "Dispatcher";
            const SchedStrategy sch =
                ist != 0 ? sn.stations[ist - 1].sched : SchedStrategy::FCFS;
            switch (sch) {
                case SchedStrategy::PS:
                case SchedStrategy::DPS:
                case SchedStrategy::GPS:
                case SchedStrategy::PSPRIO:
                case SchedStrategy::DPSPRIO:
                case SchedStrategy::GPSPRIO:
                case SchedStrategy::LPS:
                    s.server = "SharedServer";
                    break;
                case SchedStrategy::LCFSPR:
                case SchedStrategy::LCFSPRPRIO:
                case SchedStrategy::FCFSPR:
                case SchedStrategy::FCFSPRPRIO:
                case SchedStrategy::LCFSPI:
                case SchedStrategy::LCFSPIPRIO:
                case SchedStrategy::FCFSPI:
                case SchedStrategy::FCFSPIPRIO:
                case SchedStrategy::EDF:
                    s.server = "PreemptiveServer";
                    break;
                case SchedStrategy::INF:
                    s.server = "InfiniteServer";
                    break;
                case SchedStrategy::POLLING:
                    s.server = "PollingServer";
                    break;
                default:
                    s.server = "Server";
                    break;
            }
            break;
        }
    }
    return s;
}

// ---------------------------------------------------------------------------
// The metric handles, MATLAB `@@MNetwork/getAvgHandles.m`
// ---------------------------------------------------------------------------

/** The measure kinds `saveMetrics` requests, in the order it requests them. */
enum class JmtMetricKind { QLen, Util, RespT, Tput, ArvR, Tard, SysTard };

/** The JMT `type` attribute, MATLAB `MetricType.toText`. */
inline const char* jmt_metric_text(JmtMetricKind k) {
    switch (k) {
        case JmtMetricKind::QLen: return "Number of Customers";
        case JmtMetricKind::Util: return "Utilization";
        case JmtMetricKind::RespT: return "Response Time";
        case JmtMetricKind::Tput: return "Throughput";
        case JmtMetricKind::ArvR: return "Arrival Rate";
        case JmtMetricKind::Tard: return "Tardiness";
        case JmtMetricKind::SysTard: return "System Tardiness";
    }
    return "Unknown Metric";
}

/**
 * Port of the `disabled` rules in `getAvgHandles`, per (station, class).
 *
 * The rules are per kind: a Source or a Sink reports no queue length, response
 * time, utilization or tardiness but DOES report throughput and arrival rate;
 * a Fork or a Join reports no utilization; a station whose class has no service
 * process reports nothing, EXCEPT that a cache hit/miss class keeps its
 * throughput and arrival rate -- a job only ever passes through such a class,
 * it is never served in it, and disabling the measure would leave the cache's
 * hit rate unobservable.
 *
 * `has_service_tunnel` reproduces the reference's test on the SECTION object:
 * a Source, a Join and an ordinary Place have a ServiceTunnel and are therefore
 * exempt from the service-defined test entirely.
 */
template <class T>
bool jmt_metric_enabled(const qn::NetworkStruct<T>& sn, JmtMetricKind kind, std::size_t ist,
                        std::size_t r, const std::vector<bool>& is_cache_class) {
    const std::size_t ind = sn.station_to_node[ist - 1];
    const NodeType ty = sn.nodes[ind - 1].nodetype;
    const bool is_source = ty == NodeType::Source;
    const bool is_sink = ty == NodeType::Sink;
    const JmtSections<T> sec = jmt_sections(sn, ind);
    const bool tunnel = sec.server == "ServiceTunnel" || sec.server == "JobSink";
    const bool service_defined = !sn.disabled[ist - 1][r - 1];

    switch (kind) {
        case JmtMetricKind::QLen:
        case JmtMetricKind::RespT:
        case JmtMetricKind::Tard:
            if (is_source || is_sink) return false;
            return tunnel || service_defined;
        case JmtMetricKind::Util:
            if (is_source || is_sink) return false;
            if (ty == NodeType::Join || ty == NodeType::Fork) return false;
            return tunnel || service_defined;
        case JmtMetricKind::Tput:
        case JmtMetricKind::ArvR:
            if (tunnel) return true;
            return service_defined || is_cache_class[r - 1];
        case JmtMetricKind::SysTard:
            return true;
    }
    return false;
}

/** The classes a Cache switches jobs into; they keep their Tput/ArvR measures. */
template <class T>
std::vector<bool> jmt_cache_classes(const qn::NetworkStruct<T>& sn) {
    std::vector<bool> f(sn.nclasses, false);
    for (const auto& kv : sn.nodeparam) {
        const qn::CacheParam<T>& cp = kv.second;
        for (std::size_t r = 0; r < cp.hitclass.size(); ++r)
            if (cp.hitclass[r] > 0 && cp.hitclass[r] <= sn.nclasses) f[cp.hitclass[r] - 1] = true;
        for (std::size_t r = 0; r < cp.missclass.size(); ++r)
            if (cp.missclass[r] > 0 && cp.missclass[r] <= sn.nclasses) f[cp.missclass[r] - 1] = true;
    }
    return f;
}

// ---------------------------------------------------------------------------
// Small shared emitters
// ---------------------------------------------------------------------------

/** `<parameter array="true" classPath="CP" name="NAME">` */
inline xml::Element& jmt_param(xml::Element& section, const char* class_path, const char* name,
                               bool array) {
    xml::Element& p = section.add_child("parameter");
    if (array) p.set_attr("array", "true");
    p.set_attr("classPath", class_path);
    p.set_attr("name", name);
    return p;
}

/** `<parameter classPath="CP" name="NAME"><value>V</value></parameter>` */
inline void jmt_param_value(xml::Element& section, const char* class_path, const char* name,
                            const std::string& value) {
    xml::Element& p = jmt_param(section, class_path, name, false);
    p.add_text_child("value", value);
}

/** `<refClass>NAME</refClass>`, the per-class marker of an array parameter. */
inline void jmt_ref_class(xml::Element& parent, const std::string& class_name) {
    parent.add_text_child("refClass", class_name);
}

/** MATLAB `DropStrategy.toText`, the strings JMT's dropRule field expects. */
inline const char* jmt_drop_text(DropStrategy d) {
    switch (d) {
        case DropStrategy::WAITQ: return "waiting queue";
        case DropStrategy::DROP: return "drop";
        case DropStrategy::BAS: return "BAS blocking";
        case DropStrategy::BBS: return "BBS blocking";
        case DropStrategy::RSRD: return "RSRD blocking";
        case DropStrategy::RETRIAL: return "retrial";
        case DropStrategy::RETRIAL_WITH_LIMIT: return "retrial with limit";
    }
    throw InputError("SolverJMT: unrecognized drop strategy");
}

/**
 * Whether JMT's queue section can read this drop strategy at all.
 *
 * It recognizes exactly four `dropStrategies` strings -- 'drop', 'BAS blocking',
 * 'waiting queue', 'retrial' (a lookupswitch on String.hashCode in
 * `jmt/engine/NodeSections/Queue.class`; the Storage section of a Place is even
 * narrower and drops 'retrial'). An unrecognized value falls through the default arm
 * with NO flag set, so BBS, RSRD and retrial-with-limit are not approximated, they are
 * IGNORED.
 */
inline bool jmt_reads_drop(DropStrategy d) {
    return d == DropStrategy::DROP || d == DropStrategy::BAS || d == DropStrategy::WAITQ ||
           d == DropStrategy::RETRIAL;
}

/** MATLAB `HeteroSchedPolicy.toJMTText`: JMT's long descriptive identifiers. */
inline const char* jmt_hetero_text(lang::HeteroSchedPolicy p) {
    switch (p) {
        case lang::HeteroSchedPolicy::ORDER: return "Order (Assign according to order below)";
        case lang::HeteroSchedPolicy::ALIS: return "ALIS (Assign Longest Idle Server)";
        case lang::HeteroSchedPolicy::ALFS: return "ALFS (Assign Least Flexible Server)";
        case lang::HeteroSchedPolicy::FAIRNESS: return "Fairness (Move back server type when used)";
        case lang::HeteroSchedPolicy::FSF: return "FSF (Fastest Servers First)";
        case lang::HeteroSchedPolicy::RAIS: return "RAIS (Random Assignment to Idle Servers)";
    }
    return "Order (Assign according to order below)";
}

/**
 * The writer itself. One instance per exported model; it holds the derived
 * tables (connections, exportable classes, cache classes) that nearly every
 * handler reads, so they are computed once rather than per section as the
 * reference recomputes them.
 */
template <class T>
class JmtWriter {
public:
    JmtWriter(const qn::NetworkStruct<T>& sn, const JmtWriteOptions& opt)
        : sn_(sn),
          opt_(opt),
          conn_(jmt_conn_matrix(sn)),
          keep_(jmt_exportable_classes(sn)),
          cacheclass_(jmt_cache_classes(sn)) {}

    /** Port of `@@JMTIO/writeJSIM.m`; returns the serialized document. */
    std::string write_jsim() {
        std::unique_ptr<xml::Element> sim = xml::element("sim");
        save_xml_header(*sim);
        save_classes(*sim);
        for (std::size_t ind = 1; ind <= sn_.nodes.size(); ++ind) save_node(*sim, ind);
        save_metrics(*sim);
        save_links(*sim);
        save_regions(*sim);
        save_preload(*sim);
        return xml::serialize(*sim);
    }

    /**
     * The buffer-capacity refusals of `save_buffer_capacity`, as a SENTENCE
     * rather than an exception; empty when every buffer is exportable.
     *
     * ONE PREDICATE, TWO CALLERS. `save_buffer_capacity` raises it while writing
     * the JSIM document, and `jmt::jmt_method_refusal` returns it so that
     * `findSolver` never offers a jmt row on a model the writer will refuse. It
     * was reachable only from inside the writer, which is why the gate could not
     * see it and offered `jmt.jsim` on a binding buffer.
     *
     * WHAT MAKES A BUFFER BIND is not that `sn.cap` is finite: `refresh_capacity`
     * DERIVES a finite cap for every station nobody capped. It is that the cap is
     * strictly below the population that can REACH the station, and an
     * infinite-server station has no buffer at all -- the same two tests
     * `save_buffer_capacity` applies before it exports anything.
     *
     * `jmva_engine` selects the verdict, because the two engines fail a binding
     * buffer for OPPOSITE reasons. JSIM exports it whenever JMT can read the drop
     * rule, so only `assert_station_cap_exportable`'s cases go. JMVA has no
     * capacity element in its document at all -- `write_jmva` emits a station
     * type, a per-chain demand and a per-chain visit count and nothing else -- so
     * ANY binding buffer would be solved as if it were unbounded: measured on a
     * closed Delay+FCFS model, N=4, cap 2, every jmva method reported 2.19 jobs
     * at a station that can hold 2, against the exact 1.33.
     */
    std::string buffer_capacity_refusal(bool jmva_engine) const {
        for (std::size_t ist = 1; ist <= sn_.nstations; ++ist) {
            // A SOURCE AND A SINK HAVE NO BUFFER THAT CAN BIND. The Source IS the
            // external world and the Sink absorbs, so neither ever holds a job a
            // capacity could refuse, yet `refresh_capacity` writes them a row like
            // any other station. Excluded on NODE TYPE, as
            // `qn::binding_capacity_reason` excludes them, and not by name.
            const qn::NodeType ty = sn_.nodes[sn_.station_to_node[ist - 1] - 1].nodetype;
            if (ty == qn::NodeType::Source || ty == qn::NodeType::Sink) continue;
            // UNBOUNDED IS inf HERE, `Station<T>::cap` defaulting to infinity.
            // The JAR cannot: its Station.cap is an int whose "no bound" value is
            // Integer.MAX_VALUE, and refreshCapacity SUMS that sentinel across
            // the classes served, so a mixed station comes out as
            // 2147483647 + N there and needs SaveHandlers.jmtCapIsUnbounded.
            if (!std::isfinite(sn_.cap[ist - 1])) continue;
            if (sn_.cap[ist - 1] >= reachable_population(ist)) continue;
            if (!std::isfinite(sn_.stations[ist - 1].nservers)) continue;
            if (jmva_engine)
                return "SolverJMT: station '" + nname(sn_.station_to_node[ist - 1]) +
                       "' carries a finite capacity " + jmt_int(sn_.cap[ist - 1]) +
                       " that binds. The JMVA document has no capacity element at all, so the "
                       "analytical engine would solve the model as if the buffer were unbounded "
                       "and report that as the answer. Use the 'jsim' method, which exports the "
                       "buffer with its drop rule when JMT can express it, or SolverCTMC, "
                       "SolverSSA or SolverLDES";
            try {
                assert_station_cap_exportable(ist);
            } catch (const UnsupportedError& e) {
                return std::string(e.what());
            }
        }
        return std::string();
    }

private:
    const qn::NetworkStruct<T>& sn_;
    JmtWriteOptions opt_;
    std::vector<std::vector<bool>> conn_;
    std::vector<bool> keep_, cacheclass_;

    double d(const T& x) const { return num_traits<T>::to_double(x); }
    const std::string& cname(std::size_t r) const { return sn_.classes[r - 1].name; }
    const std::string& nname(std::size_t i) const { return sn_.nodes[i - 1].name; }

    /** The 1-based node indices `ind` is linked TO. */
    std::vector<std::size_t> outputs_of(std::size_t ind) const {
        std::vector<std::size_t> v;
        for (std::size_t j = 1; j <= sn_.nodes.size(); ++j)
            if (conn_[ind - 1][j - 1]) v.push_back(j);
        return v;
    }

    /** The 1-based node indices linked TO `ind`. */
    std::vector<std::size_t> inputs_of(std::size_t ind) const {
        std::vector<std::size_t> v;
        for (std::size_t j = 1; j <= sn_.nodes.size(); ++j)
            if (conn_[j - 1][ind - 1]) v.push_back(j);
        return v;
    }

    // -- header and classes ------------------------------------------------

    /** Port of `saveXMLHeader`. */
    void save_xml_header(xml::Element& sim) {
        sim.set_attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        sim.set_attr("name", opt_.file_name + ".jsimg");
        sim.set_attr("xsi:noNamespaceSchemaLocation", "SIMmodeldefinition.xsd");
        sim.set_attr("disableStatisticStop", "true");
        sim.set_attr("logDecimalSeparator", ".");
        sim.set_attr("logDelimiter", ";");
        sim.set_attr("logPath", opt_.log_path);
        sim.set_attr("logReplaceMode", "0");
        sim.set_attr("maxSamples",
                     std::isnan(opt_.max_samples) ? "10000" : jmt_int(opt_.max_samples));
        sim.set_attr("maxEvents", jmt_int(opt_.max_events));
        if (std::isfinite(opt_.max_simulated_time)) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.3f", opt_.max_simulated_time);
            sim.set_attr("maxSimulated", buf);
        }
        sim.set_attr("polling", "1.0");
        sim.set_attr("seed", jmt_int(static_cast<double>(opt_.seed)));
    }

    /**
     * Port of `saveClasses`.
     *
     * THE PRIORITY IS INVERTED. LINE orders priorities with the SMALLEST value
     * most urgent and JMT with the largest, so the exported value is
     * `max(prio) - prio(r)`. Exporting the raw number reverses the service
     * order of every priority model without any diagnostic.
     */
    void save_classes(xml::Element& sim) {
        int maxprio = 0;
        for (std::size_t r = 0; r < sn_.nclasses; ++r)
            maxprio = std::max(maxprio, sn_.classes[r].prio);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            const qn::JobClass& cl = sn_.classes[r - 1];
            xml::Element& uc = sim.add_child("userClass");
            uc.set_attr("name", cl.name);
            const bool open = !std::isfinite(cl.population);
            uc.set_attr("type", open ? "open" : "closed");
            // `sn.classdeadline` has no counterpart in this port's struct, so
            // the soft deadline is the reference's "no deadline" value. EDD and
            // EDF therefore export as the put strategies JMT does not yet
            // implement either, which is the reference's own state.
            uc.set_attr("softDeadline", "0.0");
            uc.set_attr("priority", jmt_int(static_cast<double>(maxprio - cl.prio)));
            const std::size_t refst = cl.refstat;
            const std::string refname = nname(sn_.station_to_node[refst - 1]);
            if (!open) {
                uc.set_attr("customers", jmt_int(cl.population));
                uc.set_attr("referenceSource", refname);
            } else if (sn_.disabled[refst - 1][r - 1]) {
                // An open class with no arrival process at its reference
                // station enters the model by class switching only; JMT names
                // that source 'ClassSwitch'.
                uc.set_attr("referenceSource", "ClassSwitch");
            } else {
                uc.set_attr("referenceSource", refname);
            }
        }
    }

    // -- the node/section walk ---------------------------------------------

    /** Port of the section loop of `writeJSIM`. */
    void save_node(xml::Element& sim, std::size_t ind) {
        const JmtSections<T> sec = jmt_sections(sn_, ind);
        xml::Element& node = sim.add_child("node");
        node.set_attr("name", nname(ind));
        const std::string parts[3] = {sec.input, sec.server, sec.output};
        for (int k = 0; k < 3; ++k) {
            if (parts[k].empty()) continue;
            std::string cls = parts[k];
            // `writeJSIM` promotes a Server to a PreemptiveServer for the
            // preemptive disciplines. `jmt_sections` already does that from the
            // scheduling strategy, so the promotion is not repeated here; the
            // remaining rewrites are the LINE-name to JMT-name ones.
            xml::Element& xs = node.add_child("section");
            xs.set_attr("className", cls);
            if (cls == "Buffer") {
                xs.set_attr("className", "Queue");
                save_buffer_capacity(xs, ind);
                save_drop_strategy(xs, ind);
                if (has_retrial(ind)) {
                    // With retrial JMT selects a different Queue constructor,
                    // which takes the retrial distributions BEFORE the get and
                    // put strategies and takes no impatience at all.
                    save_retrial_distributions(xs, ind);
                    save_get_strategy(xs, ind);
                    save_put_strategy(xs, ind);
                } else {
                    save_get_strategy(xs, ind);
                    save_put_strategy(xs, ind);
                    save_impatience(xs, ind);
                }
            } else if (cls == "Server") {
                save_number_of_servers(xs, ind);
                save_server_visits(xs);
                save_service_strategy(xs, ind);
                save_delay_off_strategy(xs, ind);
                // Job parallelism and heterogeneous pools. SimLoader picks the Server
                // constructor by the positional types of the parameters, so these five
                // must follow the service strategies as one block.
                save_class_parallelism(xs, ind);
                save_server_type_names(xs, ind);
                save_servers_per_type(xs, ind);
                save_server_compatibilities(xs, ind);
                save_hetero_sched_policy(xs, ind);
                warn_hetero_rates(ind);
                warn_switchover_on_non_polling(ind);
            } else if (cls == "PreemptiveServer") {
                save_number_of_servers(xs, ind);
                save_server_visits(xs);
                save_service_strategy(xs, ind);
                save_delay_off_strategy(xs, ind);
            } else if (cls == "SharedServer") {
                xs.set_attr("className", "PSServer");
                save_number_of_servers(xs, ind);
                save_server_visits(xs);
                save_service_strategy(xs, ind);
                save_delay_off_strategy(xs, ind);
                save_preemptive_strategy(xs, ind);
                save_preemptive_weights(xs, ind);
            } else if (cls == "InfiniteServer") {
                xs.set_attr("className", "Delay");
                save_service_strategy(xs, ind);
            } else if (cls == "PollingServer") {
                xs.set_attr("className", polling_server_class(ind));
                save_number_of_servers(xs, ind);
                save_server_visits(xs);
                save_service_strategy(xs, ind);
                save_switchover_strategy(xs, ind);
            } else if (cls == "RandomSource") {
                save_arrival_strategy(xs, ind);
            } else if (cls == "Dispatcher") {
                xs.set_attr("className", "Router");
                save_routing_strategy(xs, ind);
            } else if (cls == "StatelessClassSwitcher") {
                xs.set_attr("className", "ClassSwitch");
                save_class_switch_strategy(xs, ind);
            } else if (cls == "Cache") {
                save_cache_strategy(xs, ind);
            } else if (cls == "LogTunnel") {
                save_log_tunnel(xs, ind);
            } else if (cls == "Joiner") {
                xs.set_attr("className", "Join");
                save_join_strategy(xs, ind);
            } else if (cls == "Forker") {
                xs.set_attr("className", "Fork");
                save_fork_strategy(xs, ind);
            } else if (cls == "Storage") {
                save_total_capacity(xs, ind);
                save_place_capacities(xs, ind);
                save_drop_rule(xs, ind);
                save_get_strategy(xs, ind);
                save_put_strategies(xs, ind);
            } else if (cls == "Enabling") {
                save_enabling_conditions(xs, ind);
                save_inhibiting_conditions(xs, ind);
            } else if (cls == "Firing") {
                save_firing_outcomes(xs, ind);
            } else if (cls == "Timing") {
                save_mode_names(xs, ind);
                save_numbers_of_servers(xs, ind);
                save_timing_strategies(xs, ind);
                save_firing_priorities(xs, ind);
                save_firing_weights(xs, ind);
            }
            // ServiceTunnel, JobSink and Linkage carry no parameters.
        }
    }

    /** JMT's polling server class, by discipline; DECREMENTING has none. */
    const char* polling_server_class(std::size_t ind) const {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        const typename qn::NetworkStruct<T>::PollingParam pp = sn_.effective_polling(ist);
        switch (pp.ptype) {
            case PollingType::GATED: return "GatedPollingServer";
            case PollingType::EXHAUSTIVE: return "ExhaustivePollingServer";
            case PollingType::KLIMITED: return "LimitedPollingServer";
            case PollingType::DECREMENTING:
                throw UnsupportedError(
                    "SolverJMT: JMT does not support the decrementing (semiexhaustive) polling "
                    "discipline; use the LDES solver");
        }
        throw UnsupportedError("SolverJMT: unknown polling discipline");
    }

    /** True when any class of the station declares a retrial orbit. */
    bool has_retrial(std::size_t ind) const {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        if (ist == 0) return false;
        const auto it = sn_.retrialparam.find(ist);
        if (it == sn_.retrialparam.end()) return false;
        for (const lang::Distrib<T>& p : it->second.retrial_proc)
            if (!p.disabled) return true;
        return false;
    }

    /**
     * The reference WARNS that a switchover time on a non-polling queue is
     * dropped; a warning is not available here, and silently dropping a
     * declared service-order cost would answer a different model, so it is
     * refused by name.
     */
    /**
     * Port of the `writeJSIM` guard on a switchover declared away from a
     * polling server: JMT's ordinary Server has no SwitchoverStrategy, so the
     * times are WARNED ABOUT AND DROPPED, exactly as the reference does. This
     * threw instead until it was found to be the harsher rule of the two -- the
     * reference solves switchover_basic and reports a table, so a refusal here
     * left the C++ row with nothing to compare rather than with a documented
     * difference.
     */
    void warn_switchover_on_non_polling(std::size_t ind) const {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        if (ist == 0) return;
        for (const lang::Distrib<T>& s : sn_.stations[ist - 1].switchover)
            if (!s.disabled) {
                std::cerr << "[LINE] Warning: JMT does not support switchover times for "
                          << "non-polling queues. Switchover times will be ignored for node '"
                          << nname(ind) << "'." << std::endl;
                return;
            }
    }

    // -- Buffer section -----------------------------------------------------

    /**
     * Port of `saveBufferCapacity`.
     *
     * LINE's `cap` is Kendall's K, the WHOLE system capacity, and so is JMT's
     * `size`, so the number crosses unchanged. -1 is JMT's "unbounded", and it
     * is also what a capacity the population cannot REACH means: such a bound
     * can never bind, and exporting it would make JMT reject arrivals that LINE
     * admits at the instant the last job arrives.
     *
     * The test is `>=` and not `!=`: `refresh_capacity` DERIVES `sn.cap` for a
     * station the user never capped, as the sum over the classes served there
     * of the chain population, so a multi-class station gets (#classes) x N --
     * 8 on a two-class model of 4 jobs. Under `!=` only the single-class case
     * matched, and every multi-class one fell through to
     * `assert_station_cap_exportable` and was refused as a "binding" buffer
     * nobody declared. An open class carries an infinite population and so
     * makes `total` infinite, which is what keeps a DECLARED cap in a mixed
     * model refused: the open stream can fill the buffer under a closed job.
     */
    void save_buffer_capacity(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        std::string v = "-1";
        if (ist != 0 && std::isfinite(sn_.cap[ist - 1])) {
            const double total = reachable_population(ist);
            const double ns = sn_.stations[ist - 1].nservers;
            if (sn_.cap[ist - 1] < total && std::isfinite(ns)) {
                assert_station_cap_exportable(ist);
                v = jmt_int(sn_.cap[ist - 1]);
            }
        }
        jmt_param_value(section, "java.lang.Integer", "size", v);
    }

    /**
     * The most jobs that can be present at station `ist` (1-based).
     *
     * Read the way `refresh_capacity` derives the capacity itself: per CHAIN,
     * because a chain's whole population can reach a station that serves any one
     * of its classes (class switching moves jobs between them), and a chain none
     * of whose classes is served there cannot put a single job on it.
     *
     * Deliberately NOT read off the clamped per-class capacities: comparing a
     * capacity against a quantity derived from it would make every
     * user-declared buffer look non-binding. Infinite when an open chain is
     * served here, which is what the model total gave before and which sends
     * the station to `assert_station_cap_exportable`, where the open classes
     * are skipped by name.
     *
     * A class that never visits this station cannot fill it, so summing every
     * class's population made a capacity that is exactly the reachable
     * population look like a buffer -- which is what a SELF-LOOPING CLASS does.
     * See the MATLAB twin in JMTIO/saveBufferCapacity.m.
     */
    double reachable_population(std::size_t ist) const {
        double n = 0.0;
        for (std::size_t c = 0; c < sn_.nchains; ++c) {
            bool served = false;
            double chain_jobs = 0.0;
            for (std::size_t r : sn_.inchain[c]) {
                if (!sn_.disabled[ist - 1][r - 1]) served = true;
                chain_jobs += sn_.classes[r - 1].population;
            }
            if (served) n += chain_jobs;
        }
        return n;
    }

    /**
     * True when station `ist` is the RECEIVING side of a true-BAS relation for
     * class `r`, i.e. an arrival of `r` that finds `ist` full must block an
     * upstream station rather than be lost.
     *
     * LINE accepts the BAS declaration in two places -- on the blocking
     * (upstream) station, as `cqn_bas_blocking` does, or on the full
     * destination, as a model read back from JMT does -- and
     * `NetworkStruct::refresh_local_vars` resolves both into
     * `sn.isbasdestination` (BUG-83). Reading `sn.droprule` at the capped
     * station sees only the second form, which is what made SolverJMT refuse
     * the first one.
     */
    bool is_bas_destination(std::size_t ist, std::size_t r) const {
        if (ist == 0 || sn_.isbasdestination.size() < ist) return false;
        if (sn_.isbasdestination[ist - 1].size() < r) return false;
        return sn_.isbasdestination[ist - 1][r - 1];
    }

    /**
     * The JMT dropStrategy/dropRule string for station `ist`, class `r`.
     *
     * Beyond `jmt_drop_text` this resolves the two ways LINE can declare BAS
     * blocking onto the one way JMT can read it. JMT's queue section says what
     * happens to an arrival that finds THIS buffer full, so it only understands
     * the rule on the destination; a WAITQ slot that `is_bas_destination` marks
     * is therefore written out as 'BAS blocking'. A node with no buffer of its
     * own gets JMT's 'drop', as the reference does for a NaN station index.
     *
     * It also keeps the written file VALID: a strategy `jmt_reads_drop` rejects
     * is spelled 'waiting queue', JMT's own no-limit default. That substitution
     * is only ever reached where the rule cannot be consulted (infinite size, or
     * a closed capacity equal to the population): a buffer that can actually
     * fill under one of those is refused outright by
     * `assert_station_cap_exportable`.
     */
    const char* drop_strategy_text(std::size_t ist, std::size_t r) const {
        if (ist == 0) return "drop";
        const DropStrategy d = sn_.droprule[ist - 1][r - 1];
        if (d == DropStrategy::WAITQ && is_bas_destination(ist, r))
            return jmt_drop_text(DropStrategy::BAS);
        if (!jmt_reads_drop(d)) return jmt_drop_text(DropStrategy::WAITQ);
        return jmt_drop_text(d);
    }

    /**
     * Refuses a binding station capacity JMT cannot express, on two counts.
     *
     * (1) THE RULE IS ONE JMT CANNOT READ -- BBS, RSRD or retrial-with-limit;
     * see `jmt_reads_drop`. Such a value is not approximated, it is IGNORED, so
     * the capacity stops being enforced and JMT returns the unconstrained
     * answer.
     *
     * (2) THE RULE IS WAITQ AND A CLOSED CLASS CAN REACH THE LIMIT, the same
     * reason `assert_class_cap_exportable` below refuses the per-class one: JMT
     * cannot hold a blocked closed job at its upstream station. Note this is the
     * case where NO blocking rule is declared. A model that does declare BAS is
     * exported as JMT "BAS blocking", which is the same queueing model, under
     * either declaration form -- see `is_bas_destination`.
     *
     * That refusal advised expressing the limit as the STATION capacity instead,
     * and measured on 2026-08-19 the advice was wrong, and neither of the two
     * strategies a WAITQ station maps onto reproduces the UNDECLARED case:
     *
     *   waiting queue  does not enforce `size` at all. On a closed 3-queue
     *                  tandem, N=6, Exp(1) FCFS, cap 2 at Q2, JMT returned the
     *                  UNCONSTRAINED [2.03 1.99 1.98], X = 0.750, against the
     *                  exact [3.6090 0.9711 1.4199], X = 0.6522.
     *   BAS blocking   enforces it, but completes the service BEFORE blocking,
     *                  so the blocked job moves the instant room frees -- a
     *                  different queueing model, not a rounding: same fixture,
     *                  [2.871 1.373 1.756], X = 0.7126.
     *
     * With no rule declared LINE instead disables the upstream departure while
     * the destination is full, which for exponential service is repetitive
     * service (RS) and is what SolverCTMC, SolverSSA and SolverLDES all agree
     * on. So THAT model is refused rather than exported as either of the two
     * things JMT can say. See BUG-81. A declared-BAS model is a different model
     * and is exported, not refused: blocking after service is precisely what
     * JMT's "BAS blocking" does.
     *
     * That the limit CAN be reached is the caller's to establish and is not
     * retested here: `save_buffer_capacity` reaches this method only for a
     * capacity strictly below the total population, which is the one thing that
     * makes a buffer a buffer.
     */
    void assert_station_cap_exportable(std::size_t ist) const {
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (sn_.disabled[ist - 1][r - 1]) continue;  // not served here
            const DropStrategy dr = sn_.droprule[ist - 1][r - 1];
            // Unmappable for EITHER class type, so this test precedes the open-class skip
            if (!jmt_reads_drop(dr))
                throw UnsupportedError(
                    "SolverJMT: station '" + nname(sn_.station_to_node[ist - 1]) +
                    "' applies drop strategy '" + jmt_drop_text(dr) + "' to class '" + cname(r) +
                    "' and carries a finite capacity " + jmt_int(sn_.cap[ist - 1]) +
                    " it can reach. JMT's queue section reads only 'drop', 'waiting queue', "
                    "'BAS blocking' and 'retrial'; it does not approximate anything else, it "
                    "ignores it, so the capacity would stop being enforced and the run would "
                    "return the unconstrained answer. Use SolverCTMC, SolverSSA or SolverLDES");
            if (!std::isfinite(sn_.classes[r - 1].population)) continue;  // open: JMT loses it too
            if (dr != DropStrategy::WAITQ)
                continue;  // a mappable declared blocking rule is exported as itself
            if (is_bas_destination(ist, r))
                continue;  // BAS declared on the UPSTREAM station: drop_strategy_text
                           // moves it onto this one, which is where JMT reads it
            throw UnsupportedError(
                "SolverJMT: station '" + nname(sn_.station_to_node[ist - 1]) +
                "' carries a finite capacity " + jmt_int(sn_.cap[ist - 1]) +
                " that binds for the closed class '" + cname(r) +
                "'. LINE blocks a closed job that finds no room -- the upstream departure is "
                "disabled and the job stays where it is -- and no JMT drop strategy reproduces "
                "that: 'waiting queue' does not enforce the size at all, and 'BAS blocking' "
                "completes the service before blocking, which is a different queueing model. "
                "Use SolverCTMC, SolverSSA or SolverLDES, or declare DropStrategy.BAS if "
                "blocking after service is the model you want, which SolverJMT does export");
        }
    }

    /** Port of `saveDropStrategy`: the per-class rule of a full buffer. */
    void save_drop_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "java.lang.String", "dropStrategies", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            const char* txt = drop_strategy_text(ist, r);
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.String");
            sp.set_attr("name", "dropStrategy");
            sp.add_text_child("value", txt);
        }
    }

    /**
     * Port of `saveGetStrategy`: FCFS everywhere except a polling queue, whose
     * discipline decides which of JMT's three polling get strategies serves it.
     */
    void save_get_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        const bool polling = sn_.nodes[ind - 1].nodetype == NodeType::Queue && ist != 0 &&
                             sn_.stations[ist - 1].sched == SchedStrategy::POLLING;
        if (!polling) {
            xml::Element& p = section.add_child("parameter");
            p.set_attr("classPath", "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy");
            p.set_attr("name", "FCFSstrategy");
            return;
        }
        const typename qn::NetworkStruct<T>::PollingParam pp = sn_.effective_polling(ist);
        xml::Element& p = section.add_child("parameter");
        switch (pp.ptype) {
            case PollingType::GATED:
                p.set_attr("classPath",
                           "jmt.engine.NetStrategies.QueueGetStrategies.GatedPollingGetStrategy");
                break;
            case PollingType::EXHAUSTIVE:
                p.set_attr(
                    "classPath",
                    "jmt.engine.NetStrategies.QueueGetStrategies.ExhaustivePollingGetStrategy");
                break;
            case PollingType::KLIMITED: {
                p.set_attr("classPath",
                           "jmt.engine.NetStrategies.QueueGetStrategies.LimitedPollingGetStrategy");
                xml::Element& k = p.add_child("subParameter");
                k.set_attr("classPath", "java.lang.Integer");
                k.set_attr("name", "pollingKValue");
                k.add_text_child("value", jmt_int(static_cast<double>(pp.pk)));
                break;
            }
            case PollingType::DECREMENTING:
                throw UnsupportedError(
                    "SolverJMT: JMT does not support the decrementing (semiexhaustive) polling "
                    "discipline; use the LDES solver");
        }
        p.set_attr("name", "FCFSstrategy");
    }

    /**
     * Port of `savePutStrategy`: the discipline expressed as WHERE an arrival
     * is inserted in the buffer.
     *
     * JMT has no scheduling-strategy field: FCFS is a tail insertion, LCFS a
     * head insertion, SJF/SRPT an ordered one, and a preemptive discipline its
     * own put strategy. Everything unlisted -- PS above all, whose sharing is
     * expressed by the PSServer section instead -- is a tail insertion.
     */
    void save_put_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p =
            jmt_param(section, "jmt.engine.NetStrategies.QueuePutStrategy", "QueuePutStrategy", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            const char* nm = "TailStrategy";
            if (ist != 0) {
                switch (sn_.stations[ist - 1].sched) {
                    case SchedStrategy::SIRO: nm = "RandStrategy"; break;
                    case SchedStrategy::LJF: nm = "LJFStrategy"; break;
                    case SchedStrategy::SJF: nm = "SJFStrategy"; break;
                    case SchedStrategy::LEPT: nm = "LEPTStrategy"; break;
                    case SchedStrategy::SEPT: nm = "SEPTStrategy"; break;
                    case SchedStrategy::LCFS: nm = "HeadStrategy"; break;
                    case SchedStrategy::LCFSPRIO: nm = "HeadStrategyPriority"; break;
                    case SchedStrategy::LCFSPR: nm = "LCFSPRStrategy"; break;
                    case SchedStrategy::LCFSPI: nm = "LCFSPIStrategy"; break;
                    case SchedStrategy::LCFSPRPRIO: nm = "LCFSPRStrategyPriority"; break;
                    case SchedStrategy::LCFSPIPRIO: nm = "LCFSPIStrategyPriority"; break;
                    case SchedStrategy::FCFSPR: nm = "FCFSPRStrategy"; break;
                    case SchedStrategy::FCFSPI: nm = "FCFSPIStrategy"; break;
                    case SchedStrategy::FCFSPRPRIO: nm = "FCFSPRStrategyPriority"; break;
                    case SchedStrategy::FCFSPIPRIO: nm = "FCFSPIStrategyPriority"; break;
                    case SchedStrategy::HOL: nm = "TailStrategyPriority"; break;
                    case SchedStrategy::EDD: nm = "EDDStrategy"; break;
                    case SchedStrategy::EDF: nm = "EDFStrategy"; break;
                    case SchedStrategy::SRPT: nm = "SRPTStrategy"; break;
                    case SchedStrategy::SRPTPRIO: nm = "SRPTStrategyPriority"; break;
                    default: nm = "TailStrategy"; break;
                }
            }
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath",
                        std::string("jmt.engine.NetStrategies.QueuePutStrategies.") + nm);
            sp.set_attr("name", nm);
        }
    }

    // -- Server section -----------------------------------------------------

    /**
     * Port of `saveNumberOfServers`.
     *
     * LPS EXPORTS AS ONE SERVER. Its admission limit is not a server count but
     * a cap on the number in service, which `save_regions` expresses as an
     * implicit finite capacity region; exporting the limit as the server count
     * would give a c-server FCFS queue instead of limited processor sharing.
     * A load-dependent station exports `max(nservers, max lldscaling)`, since
     * `min(1:N, c)` scaling IS a c-server queue.
     */
    void save_number_of_servers(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        double maxjobs = 1.0;
        if (ist != 0) {
            if (sn_.stations[ist - 1].sched == SchedStrategy::LPS) {
                maxjobs = 1.0;
            } else {
                maxjobs = sn_.stations[ist - 1].nservers;
                for (const T& s : sn_.stations[ist - 1].lldscaling)
                    maxjobs = std::max(maxjobs, d(s));
            }
        }
        jmt_param_value(section, "java.lang.Integer", "maxJobs", jmt_int(maxjobs));
    }

    /** Port of `saveServerVisits`: one visit per class, always. */
    void save_server_visits(xml::Element& section) {
        xml::Element& p = jmt_param(section, "java.lang.Integer", "numberOfVisits", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Integer");
            sp.set_attr("name", "numberOfVisits");
            sp.add_text_child("value", "1");
        }
    }

    /**
     * Server pools and job parallelism of a node, as the JMT Server section needs them.
     *
     * JMT's Server section takes classParallelism, serverNames, serversPerServerType,
     * serverCompatibilities and schedulingPolicy as one positional block of its
     * constructor (`jmt.engine.NodeSections.Server`), so the five are emitted
     * together or not at all, and always after the service strategies. A station
     * declaring parallelism alone is therefore given one synthetic pool holding all
     * of its servers, since the pool counts, not numberOfServers, size the server
     * pool once any pool is declared.
     */
    struct ServerPools {
        bool present = false;
        std::vector<std::string> names;
        std::vector<double> counts;
        std::vector<std::vector<bool>> compat;  ///< [type][class]
        lang::HeteroSchedPolicy policy = lang::HeteroSchedPolicy::ORDER;
        std::vector<std::size_t> parallelism;   ///< per class, one-based class order
    };

    ServerPools server_pools(std::size_t ind) const {
        ServerPools pools;
        const std::size_t ist = sn_.nodes[ind - 1].station;
        if (ist == 0) return pools;
        const auto& st = sn_.stations[ist - 1];
        const bool has_types = !st.server_types.empty();
        bool has_parallelism = false;
        for (std::size_t n : st.server_parallelism) {
            if (n > 1) { has_parallelism = true; break; }
        }
        if (!has_types && !has_parallelism) return pools;

        pools.present = true;
        pools.parallelism.assign(sn_.nclasses, 1);
        for (std::size_t r = 0; r < sn_.nclasses && r < st.server_parallelism.size(); ++r) {
            pools.parallelism[r] = st.server_parallelism[r] < 1 ? 1 : st.server_parallelism[r];
        }
        if (has_types) {
            pools.policy = st.hetero_policy;
            for (const auto& pool : st.server_types) {
                pools.names.push_back(pool.name);
                pools.counts.push_back(pool.count);
                std::vector<bool> row(sn_.nclasses, true);
                for (std::size_t r = 0; r < sn_.nclasses; ++r) {
                    // An EMPTY compatibility row means "every class", the constructor
                    // default; JMT has no such shorthand, so it is expanded here.
                    row[r] = pool.compatible.empty()
                                 ? true
                                 : (r < pool.compatible.size() ? pool.compatible[r] : false);
                }
                pools.compat.push_back(row);
            }
        } else {
            pools.names.push_back(sn_.nodes[ind - 1].name + " - Server Type 1");
            pools.counts.push_back(st.nservers);
            pools.compat.push_back(std::vector<bool>(sn_.nclasses, true));
        }
        return pools;
    }

    /** Port of `saveClassParallelism`: `Server.serverNumRequired`, per class. */
    void save_class_parallelism(xml::Element& section, std::size_t ind) {
        const ServerPools pools = server_pools(ind);
        if (!pools.present) return;
        xml::Element& p = jmt_param(section, "java.lang.Integer", "classParallelism", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Integer");
            sp.set_attr("name", "serverParallelism");
            sp.add_text_child("value", jmt_int(static_cast<double>(pools.parallelism[r - 1])));
        }
    }

    /** Port of `saveServerTypeNames`. */
    void save_server_type_names(xml::Element& section, std::size_t ind) {
        const ServerPools pools = server_pools(ind);
        if (!pools.present) return;
        xml::Element& p = jmt_param(section, "java.lang.String", "serverNames", true);
        for (const std::string& name : pools.names) {
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.String");
            sp.set_attr("name", "serverTypesNames");
            sp.add_text_child("value", name);
        }
    }

    /** Port of `saveServersPerType`. */
    void save_servers_per_type(xml::Element& section, std::size_t ind) {
        const ServerPools pools = server_pools(ind);
        if (!pools.present) return;
        xml::Element& p = jmt_param(section, "java.lang.Integer", "serversPerServerType", true);
        for (double count : pools.counts) {
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Integer");
            sp.set_attr("name", "serverTypesNumOfServers");
            sp.add_text_child("value", jmt_int(count));
        }
    }

    /** Port of `saveServerCompatibilities`. */
    void save_server_compatibilities(xml::Element& section, std::size_t ind) {
        const ServerPools pools = server_pools(ind);
        if (!pools.present) return;
        xml::Element& p = jmt_param(section, "java.lang.Object", "serverCompatibilities", true);
        for (const std::vector<bool>& row : pools.compat) {
            xml::Element& tn = p.add_child("subParameter");
            tn.set_attr("array", "true");
            tn.set_attr("classPath", "java.lang.Boolean");
            tn.set_attr("name", "serverTypesCompatibilities");
            for (std::size_t r = 0; r < sn_.nclasses; ++r) {
                xml::Element& cn = tn.add_child("subParameter");
                cn.set_attr("classPath", "java.lang.Boolean");
                cn.set_attr("name", "compatibilities");
                cn.add_text_child("value", row[r] ? "true" : "false");
            }
        }
    }

    /** Port of `saveHeteroSchedPolicy`, via `HeteroSchedPolicy.toJMTText`. */
    void save_hetero_sched_policy(xml::Element& section, std::size_t ind) {
        const ServerPools pools = server_pools(ind);
        if (!pools.present) return;
        jmt_param_value(section, "java.lang.String", "schedulingPolicy",
                        jmt_hetero_text(pools.policy));
    }

    /**
     * Warns that per-server-type service rates cannot reach the JMT engine.
     *
     * JMT keys the ServiceStrategy array of a station by refClass, so its loader
     * (`jmt.engine.simEngine.SimLoader`) keeps one strategy per class however many
     * (type, class) entries are written, and every pool of a station ends up
     * serving at the class rate. Pool sizes, class compatibilities and the
     * assignment policy do cross; per-pool service laws do not.
     */
    void warn_hetero_rates(std::size_t ind) const {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        if (ist == 0) return;
        const auto& pools = sn_.stations[ist - 1].server_types;
        if (pools.size() < 2) return;
        bool first_set = false, distinct = false;
        double first = 0.0;
        for (const auto& pool : pools) {
            for (const lang::Distrib<T>& sd : pool.service) {
                if (sd.disabled) continue;
                const double mean = num_traits<T>::to_double(sd.mean);
                if (!(mean > 0)) continue;
                if (!first_set) { first = mean; first_set = true; }
                else if (std::abs(mean - first) > 1e-12) { distinct = true; }
            }
        }
        if (!distinct) return;
        std::cerr << "[LINE] Warning: JMT keys service strategies by job class, so the "
                  << "per-server-type service rates of station '" << sn_.nodes[ind - 1].name
                  << "' cannot be exported; every pool will serve at the class service rate. "
                  << "Use the LDES or CTMC solver for per-type rates." << std::endl;
    }

    /**
     * Port of `saveServiceStrategy`.
     *
     * A pair the class never visits becomes a DisabledServiceTimeStrategy and
     * an Immediate becomes a ZeroServiceTimeStrategy; neither carries a
     * distribution. Everything else goes through the shared emitter.
     */
    void save_service_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "jmt.engine.NetStrategies.ServiceStrategy",
                                    "ServiceStrategy", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& sts = p.add_child("subParameter");
            const JmtDistView<T> v =
                ist == 0 ? JmtDistView<T>() : jmt_dist_view(sn_.service[ist - 1][r - 1]);
            if (v.type == lang::ProcessType::DISABLED) {
                sts.set_attr("classPath",
                             "jmt.engine.NetStrategies.ServiceStrategies."
                             "DisabledServiceTimeStrategy");
                sts.set_attr("name", "DisabledServiceTimeStrategy");
                continue;
            }
            if (v.type == lang::ProcessType::IMMEDIATE) {
                sts.set_attr("classPath",
                             "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                sts.set_attr("name", "ZeroServiceTimeStrategy");
                continue;
            }
            sts.set_attr("classPath",
                         "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            sts.set_attr("name", "ServiceTimeStrategy");
            jmt_append_distribution(sts, v, "SolverJMT (service)");
        }
    }

    /**
     * Port of `saveArrivalStrategy`.
     *
     * A CLOSED class has no arrival process at the Source and is exported as a
     * ServiceTimeStrategy with a literal `null` body, which is how JMT spells
     * "this class does not arrive here". So is an open class the Source
     * disables -- one that enters by class switching only.
     */
    void save_arrival_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "jmt.engine.NetStrategies.ServiceStrategy",
                                    "ServiceStrategy", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& sts = p.add_child("subParameter");
            sts.set_attr("classPath",
                         "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            sts.set_attr("name", "ServiceTimeStrategy");
            const bool closed = std::isfinite(sn_.classes[r - 1].population);
            const JmtDistView<T> v =
                ist == 0 ? JmtDistView<T>() : jmt_dist_view(sn_.service[ist - 1][r - 1]);
            if (closed || v.type == lang::ProcessType::DISABLED) {
                sts.add_text_child("value", "null");
                continue;
            }
            if (v.type == lang::ProcessType::IMMEDIATE) {
                sts.set_attr("classPath",
                             "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                sts.set_attr("name", "ZeroServiceTimeStrategy");
                continue;
            }
            jmt_append_distribution(sts, v, "SolverJMT (arrival)");
        }
    }

    /** Port of `savePreemptiveStrategy`: which PS variant a PSServer runs. */
    void save_preemptive_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "jmt.engine.NetStrategies.PSStrategy", "PSStrategy",
                                    true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            const char* nm = nullptr;
            if (ist != 0) switch (sn_.stations[ist - 1].sched) {
                    // LPS shares the server evenly among those admitted, so its
                    // in-service discipline IS EPS; the limit is the region.
                    case SchedStrategy::PS:
                    case SchedStrategy::LPS: nm = "EPSStrategy"; break;
                    case SchedStrategy::DPS: nm = "DPSStrategy"; break;
                    case SchedStrategy::GPS: nm = "GPSStrategy"; break;
                    case SchedStrategy::PSPRIO: nm = "EPSStrategyPriority"; break;
                    case SchedStrategy::DPSPRIO: nm = "DPSStrategyPriority"; break;
                    case SchedStrategy::GPSPRIO: nm = "GPSStrategyPriority"; break;
                    default: nm = nullptr; break;
                }
            xml::Element& sp = p.add_child("subParameter");
            if (nm != nullptr) {
                sp.set_attr("classPath", std::string("jmt.engine.NetStrategies.PSStrategies.") + nm);
                sp.set_attr("name", nm);
            }
        }
    }

    /** Port of `savePreemptiveWeights`: the DPS/GPS share of each class. */
    void save_preemptive_weights(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "java.lang.Double", "serviceWeights", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            double w = 0.0;
            if (ist != 0 && sn_.stations[ist - 1].schedparam.size() >= r)
                w = d(sn_.stations[ist - 1].schedparam[r - 1]);
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Double");
            sp.set_attr("name", "serviceWeight");
            sp.add_text_child("value", jmt_num(w));
        }
    }

    /**
     * Port of `saveDelayOffStrategy`: the setup and delay-off times of a server
     * that powers down when idle.
     *
     * JMT casts each per-class entry to `ServiceStrategy[]` and reads element
     * [0], so the strategy sits inside a single-element array rather than
     * directly under the class. A class with no time declared gets a
     * deterministic zero, not an omission: an absent entry would leave JMT's
     * array short and shift every later class's setup time onto the wrong one.
     */
    void save_delay_off_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        if (ist == 0) return;
        const auto it = sn_.setupparam.find(ist);
        if (it == sn_.setupparam.end()) return;
        const qn::SetupDelayOffParam<T>& sp = it->second;
        bool any = false;
        for (const lang::Distrib<T>& s : sp.setup)
            if (!s.disabled) any = true;
        if (!any) return;

        const char* names[2] = {"delayOffTime", "setUpTime"};
        for (int which = 0; which < 2; ++which) {
            const std::vector<lang::Distrib<T>>& tab = which == 0 ? sp.delayoff : sp.setup;
            xml::Element& p = jmt_param(section, "java.lang.Object", names[which], true);
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                jmt_ref_class(p, cname(r));
                xml::Element& row = p.add_child("subParameter");
                row.set_attr("array", "true");
                row.set_attr("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
                row.set_attr("name", names[which]);
                xml::Element& sts = row.add_child("subParameter");
                sts.set_attr("classPath",
                             "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                sts.set_attr("name", "ServiceTimeStrategy");
                if (r <= tab.size() && !tab[r - 1].disabled)
                    append_simple_distribution(sts, tab[r - 1]);
                else
                    append_zero_time(sts);
            }
        }
    }

    /** `appendZeroTimeXml`: a deterministic zero. */
    void append_zero_time(xml::Element& parent) {
        xml::Element& dn = parent.add_child("subParameter");
        dn.set_attr("classPath", "jmt.engine.random.DeterministicDistr");
        dn.set_attr("name", "Deterministic");
        xml::Element& par = parent.add_child("subParameter");
        par.set_attr("classPath", "jmt.engine.random.DeterministicDistrPar");
        par.set_attr("name", "distrPar");
        jmt_scalar(par, "java.lang.Double", "t", "0.0");
    }

    /**
     * `appendDistributionXml`: the REDUCED emitter the setup/delay-off pair
     * uses -- Exp, Erlang, Det, Immediate, and everything else as a
     * deterministic time equal to its mean.
     *
     * The fallback is the reference's and is kept: a setup time is a small
     * fixed overhead in every model that has one, so replacing an exotic law by
     * its mean there is a documented simplification rather than a silent one,
     * and refusing would reject models JMT can otherwise run.
     */
    void append_simple_distribution(xml::Element& parent, const lang::Distrib<T>& dist) {
        if (dist.type == lang::ProcessType::IMMEDIATE) {
            append_zero_time(parent);
            return;
        }
        const double mean = d(dist.mean);
        if (dist.type == lang::ProcessType::EXP) {
            xml::Element& dn = parent.add_child("subParameter");
            dn.set_attr("classPath", "jmt.engine.random.Exponential");
            dn.set_attr("name", "Exponential");
            xml::Element& par = parent.add_child("subParameter");
            par.set_attr("classPath", "jmt.engine.random.ExponentialPar");
            par.set_attr("name", "distrPar");
            jmt_double(par, "lambda", mean > 0.0 ? 1.0 / mean : 0.0);
            return;
        }
        if (dist.type == lang::ProcessType::ERLANG) {
            const std::size_t ph = dist.phases();
            xml::Element& dn = parent.add_child("subParameter");
            dn.set_attr("classPath", "jmt.engine.random.Erlang");
            dn.set_attr("name", "Erlang");
            xml::Element& par = parent.add_child("subParameter");
            par.set_attr("classPath", "jmt.engine.random.ErlangPar");
            par.set_attr("name", "distrPar");
            jmt_double(par, "alpha", mean > 0.0 ? static_cast<double>(ph) / mean : 0.0);
            jmt_scalar(par, "java.lang.Long", "r", jmt_int(static_cast<double>(ph)));
            return;
        }
        xml::Element& dn = parent.add_child("subParameter");
        dn.set_attr("classPath", "jmt.engine.random.DeterministicDistr");
        dn.set_attr("name", "Deterministic");
        xml::Element& par = parent.add_child("subParameter");
        par.set_attr("classPath", "jmt.engine.random.DeterministicDistrPar");
        par.set_attr("name", "distrPar");
        jmt_double(par, "t", mean);
    }

    // -- Dispatcher / ClassSwitch / Fork / Join ------------------------------

    /**
     * Port of `saveRoutingStrategy`.
     *
     * A CLASS SWITCH NODE IS EXPORTED AS RANDOM ROUTING. It has exactly one
     * outgoing link, so uniform routing over it is the same routing, and it
     * avoids reading `rt` at a node whose class-switch mass the refresh has
     * already folded into the edges around it.
     */
    void save_routing_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t K = sn_.nclasses;
        const bool is_cs = sn_.nodes[ind - 1].nodetype == NodeType::ClassSwitch;
        xml::Element& p =
            jmt_param(section, "jmt.engine.NetStrategies.RoutingStrategy", "RoutingStrategy", true);
        for (std::size_t r = 1; r <= K; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            const RoutingStrategy rs =
                is_cs ? RoutingStrategy::RAND
                      : (sn_.nodes[ind - 1].routing.size() >= r ? sn_.nodes[ind - 1].routing[r - 1]
                                                                : RoutingStrategy::PROB);
            xml::Element& sp = p.add_child("subParameter");
            switch (rs) {
                case RoutingStrategy::RAND:
                    sp.set_attr("classPath",
                                "jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy");
                    sp.set_attr("name", "Random");
                    break;
                case RoutingStrategy::RROBIN:
                    sp.set_attr("classPath",
                                "jmt.engine.NetStrategies.RoutingStrategies.RoundRobinStrategy");
                    sp.set_attr("name", "Round Robin");
                    break;
                case RoutingStrategy::JSQ:
                    sp.set_attr("classPath",
                                "jmt.engine.NetStrategies.RoutingStrategies."
                                "ShortestQueueLengthRoutingStrategy");
                    sp.set_attr("name", "Join the Shortest Queue (JSQ)");
                    break;
                case RoutingStrategy::SQ: {
                    sp.set_attr("classPath",
                                "jmt.engine.NetStrategies.RoutingStrategies.PowerOfKRoutingStrategy");
                    sp.set_attr("name", "Power of k");
                    const int dpar = sn_.nodes[ind - 1].routing_param.size() >= r
                                         ? sn_.nodes[ind - 1].routing_param[r - 1]
                                         : 0;
                    jmt_scalar(sp, "java.lang.Integer", "k", jmt_int(static_cast<double>(dpar)));
                    // Always false: LINE implements SQ(d) only. JMT's
                    // withMemory selects Anselmi & Dufour SQ(d,N), a different
                    // policy.
                    jmt_scalar(sp, "java.lang.Boolean", "withMemory", "false");
                    break;
                }
                case RoutingStrategy::WRROBIN: {
                    sp.set_attr("classPath",
                                "jmt.engine.NetStrategies.RoutingStrategies."
                                "WeightedRoundRobinStrategy");
                    sp.set_attr("name", "Weighted Round Robin");
                    xml::Element& arr =
                        jmt_param_sub(sp, "jmt.engine.NetStrategies.RoutingStrategies.WeightEntry",
                                      "WeightEntryArray");
                    const std::map<std::size_t, double>& w =
                        sn_.nodes[ind - 1].routing_weights.size() >= r
                            ? sn_.nodes[ind - 1].routing_weights[r - 1]
                            : empty_weights_;
                    for (std::size_t j : outputs_of(ind)) {
                        const auto wi = w.find(j);
                        xml::Element& e = arr.add_child("subParameter");
                        e.set_attr("classPath",
                                   "jmt.engine.NetStrategies.RoutingStrategies.WeightEntry");
                        e.set_attr("name", "WeightEntry");
                        jmt_scalar(e, "java.lang.String", "stationName", nname(j));
                        jmt_scalar(e, "java.lang.Integer", "weight",
                                   jmt_int(wi == w.end() ? 0.0 : wi->second));
                    }
                    break;
                }
                case RoutingStrategy::PROB: {
                    sp.set_attr("classPath",
                                "jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy");
                    sp.set_attr("name", "Probabilities");
                    xml::Element& arr = jmt_param_sub(sp, "jmt.engine.random.EmpiricalEntry",
                                                      "EmpiricalEntryArray");
                    for (std::size_t j : outputs_of(ind)) {
                        const double pr = route_node(ind, r, j, r);
                        if (!(pr > 0.0)) continue;
                        xml::Element& e = arr.add_child("subParameter");
                        e.set_attr("classPath", "jmt.engine.random.EmpiricalEntry");
                        e.set_attr("name", "EmpiricalEntry");
                        jmt_scalar(e, "java.lang.String", "stationName", nname(j));
                        jmt_scalar(e, "java.lang.Double", "probability", jmt_fmt(pr));
                    }
                    break;
                }
                default:
                    sp.set_attr(
                        "classPath",
                        "jmt.engine.NetStrategies.RoutingStrategies.DisabledRoutingStrategy");
                    sp.set_attr("name", "Random");
                    break;
            }
        }
    }

    /** `<subParameter array="true" classPath=CP name=NAME>` */
    static xml::Element& jmt_param_sub(xml::Element& parent, const char* cp, const char* name) {
        xml::Element& e = parent.add_child("subParameter");
        e.set_attr("array", "true");
        e.set_attr("classPath", cp);
        e.set_attr("name", name);
        return e;
    }

    /** `sn.rtnodes((i-1)K+r, (j-1)K+s)`, guarded against an unbuilt matrix. */
    double route_node(std::size_t i, std::size_t r, std::size_t j, std::size_t s) const {
        const std::size_t K = sn_.nclasses;
        const std::size_t a = (i - 1) * K + (r - 1), b = (j - 1) * K + (s - 1);
        if (sn_.rtnodes.rows() <= a || sn_.rtnodes.cols() <= b) return 0.0;
        return d(sn_.rtnodes(a, b));
    }

    /**
     * Port of `saveClassSwitchStrategy`: the (K x K) switching matrix as JMT
     * reads it, each cell the mass class r leaves this node with as class s.
     *
     * The cell is the SUM over the node's outgoing links of `rtnodes`, not the
     * declared `csmatrix`: the refresh has already multiplied the switch by the
     * routing, and the sum recovers the switching probability independently of
     * how the mass was split over the links.
     */
    void save_class_switch_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t K = sn_.nclasses;
        const std::vector<std::size_t> jset = outputs_of(ind);
        xml::Element& p = jmt_param(section, "java.lang.Object", "matrix", true);
        for (std::size_t r = 1; r <= K; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& row = p.add_child("subParameter");
            row.set_attr("array", "true");
            row.set_attr("classPath", "java.lang.Float");
            row.set_attr("name", "row");
            for (std::size_t s = 1; s <= K; ++s) {
                if (!keep_[s - 1]) continue;
                jmt_ref_class(row, cname(s));
                double acc = 0.0;
                for (std::size_t j : jset) acc += route_node(ind, r, j, s);
                xml::Element& cell = row.add_child("subParameter");
                cell.set_attr("classPath", "java.lang.Float");
                cell.set_attr("name", "cell");
                cell.add_text_child("value", jmt_fmt(acc));
            }
        }
    }

    /**
     * Port of `saveForkStrategy`.
     *
     * ONE OutPath ENTRY PER BRANCH. The reference used to build the entry
     * inside a loop over the connected nodes but append it outside, so only the
     * LAST link survived, and the JAR and this port transcribed the same shape.
     * It was harmless only because `isSimplifiedFork` is true, which makes JMT
     * send one job down every outgoing link and ignore the branch list -- and it
     * stops being harmless the moment a branch carries its own probability or
     * its own jobs-per-link. All four codebases now emit the full list.
     */
    void save_fork_strategy(xml::Element& section, std::size_t ind) {
        const qn::ForkParam<T>* fk = sn_.fork_param_of(ind);
        const double fan_out = sn_.nodes[ind - 1].tasks_per_link;
        jmt_param_value(section, "java.lang.Integer", "jobsPerLink", jmt_int(fan_out));
        jmt_param_value(section, "java.lang.Integer", "block", "-1");

        // isSimplifiedFork lets JMT ignore the branch list and send one job down
        // every link. That is only the same model when every branch is certain
        // and carries the same number of tasks, so a variable forking level
        // switches it off and makes JMT read the per-branch entries below.
        bool simplified = true;
        if (fk != 0) {
            for (std::size_t k = 0; k < fk->fan_out_link.rows() && simplified; ++k)
                for (std::size_t r = 0; r < fk->fan_out_link.cols() && simplified; ++r) {
                    const double p = num_traits<T>::to_double(fk->fan_out_prob(k, r));
                    if (p == 0.0) continue;  // link not taken
                    if (p != 1.0) simplified = false;
                    if (num_traits<T>::to_double(fk->fan_out_link(k, r)) != fan_out)
                        simplified = false;
                    if (!fk->fan_out_dist[k][r].disabled) simplified = false;
                }
        }
        jmt_param_value(section, "java.lang.Boolean", "isSimplifiedFork",
                        simplified ? "true" : "false");

        xml::Element& p =
            jmt_param(section, "jmt.engine.NetStrategies.ForkStrategy", "ForkStrategy", true);
        const std::vector<std::size_t> outs = outputs_of(ind);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& cs = p.add_child("subParameter");
            cs.set_attr("classPath", "jmt.engine.NetStrategies.ForkStrategies.ProbabilitiesFork");
            cs.set_attr("name", "Branch Probabilities");
            xml::Element& arr = jmt_param_sub(
                cs, "jmt.engine.NetStrategies.ForkStrategies.OutPath", "EmpiricalEntryArray");
            const RoutingStrategy rs = sn_.nodes[ind - 1].routing.size() >= r
                                           ? sn_.nodes[ind - 1].routing[r - 1]
                                           : RoutingStrategy::PROB;
            if (rs != RoutingStrategy::PROB || outs.empty()) continue;
            for (std::size_t oi = 0; oi < outs.size(); ++oi) {
                xml::Element& entry = arr.add_child("subParameter");
                entry.set_attr("classPath", "jmt.engine.NetStrategies.ForkStrategies.OutPath");
                entry.set_attr("name", "OutPathEntry");
                xml::Element& unit = entry.add_child("subParameter");
                unit.set_attr("classPath", "jmt.engine.random.EmpiricalEntry");
                unit.set_attr("name", "outUnitProbability");
                jmt_scalar(unit, "java.lang.String", "stationName", nname(outs[oi]));

                // branch activation probability: JMT's outUnitProbability
                const std::size_t k0 = outs[oi] - 1, r0 = r - 1;
                const bool has_fan = fk != 0;
                const double branch_p =
                    has_fan ? num_traits<T>::to_double(fk->fan_out_prob(k0, r0)) : 1.0;
                jmt_scalar(unit, "java.lang.Double", "probability", jmt_fmt(branch_p));

                // JobsPerLinkDis is an EmpiricalEntry ARRAY: one entry per point
                // of the jobs-per-link distribution. A deterministic fork emits
                // the single degenerate entry it always did.
                std::vector<double> pts, prs;
                if (has_fan && !fk->fan_out_dist[k0][r0].disabled) {
                    const lang::Distrib<T>& d = fk->fan_out_dist[k0][r0];
                    double tot = 0.0;
                    for (std::size_t e = 0; e < d.params.size(); ++e)
                        tot += num_traits<T>::to_double(d.params[e]);
                    for (std::size_t e = 0; e < d.params.size(); ++e) {
                        pts.push_back(d.trace.empty()
                                          ? static_cast<double>(e + 1)
                                          : num_traits<T>::to_double(d.trace[e]));
                        prs.push_back(num_traits<T>::to_double(d.params[e]) / tot);
                    }
                } else if (has_fan) {
                    pts.push_back(num_traits<T>::to_double(fk->fan_out_link(k0, r0)));
                    prs.push_back(1.0);
                } else {
                    pts.push_back(fan_out);
                    prs.push_back(1.0);
                }

                xml::Element& jpl =
                    jmt_param_sub(entry, "jmt.engine.random.EmpiricalEntry", "JobsPerLinkDis");
                for (std::size_t e = 0; e < pts.size(); ++e) {
                    xml::Element& jple = jpl.add_child("subParameter");
                    jple.set_attr("classPath", "jmt.engine.random.EmpiricalEntry");
                    jple.set_attr("name", "EmpiricalEntry");
                    jmt_scalar(jple, "java.lang.String", "numbers", jmt_int(pts[e]));
                    jmt_scalar(jple, "java.lang.Double", "probability", jmt_fmt(prs[e]));
                }
            }
        }
    }

    /**
     * Port of `saveJoinStrategy`.
     *
     * `numRequired` is -1 for a standard join -- JMT's "every sibling" -- and
     * the quorum for a partial one. The reference reads two separate fields
     * (`fanIn` and `joinRequired`); this port's `JoinDecl` carries the quorum
     * once, with 0 meaning "every sibling", which is the same information.
     */
    void save_join_strategy(xml::Element& section, std::size_t ind) {
        const auto it = sn_.joindecl.find(ind);
        const bool partial =
            it != sn_.joindecl.end() && it->second.strategy == JoinStrategy::PARTIAL;
        const double quorum = it != sn_.joindecl.end() ? it->second.quorum : 0.0;
        xml::Element& p =
            jmt_param(section, "jmt.engine.NetStrategies.JoinStrategy", "JoinStrategy", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            jmt_ref_class(p, cname(r));
            xml::Element& sp = p.add_child("subParameter");
            if (partial) {
                sp.set_attr("classPath", "jmt.engine.NetStrategies.JoinStrategies.PartialJoin");
                sp.set_attr("name", "Quorum");
            } else {
                sp.set_attr("classPath", "jmt.engine.NetStrategies.JoinStrategies.NormalJoin");
                sp.set_attr("name", "Standard Join");
            }
            const double req = partial ? quorum : (quorum > 0.0 ? quorum : -1.0);
            jmt_scalar(sp, "java.lang.Integer", "numRequired", jmt_int(req));
        }
    }

    /** Port of `saveLogTunnel`: the ten fields JMT's LogTunnel section takes. */
    void save_log_tunnel(xml::Element& section, std::size_t ind) {
        const qn::NodeDef::LoggerParam& lg = sn_.nodes[ind - 1].logger;
        std::string path = lg.file_path.empty() ? sn_.log_path : lg.file_path;
        if (!path.empty() && path[path.size() - 1] != '/') path += '/';
        const char* bools[7] = {"logExecTimestamp", "logLoggerName", "logTimeStamp",
                                "logJobID",        "logJobClass",   "logTimeSameClass",
                                "logTimeAnyClass"};
        const bool vals[7] = {lg.start_time, lg.logger_name,     lg.timestamp,      lg.job_id,
                              lg.job_class,  lg.time_same_class, lg.time_any_class};
        jmt_param_value(section, "java.lang.String", "logfileName", lg.file_name);
        jmt_param_value(section, "java.lang.String", "logfilePath", path);
        for (int j = 0; j < 7; ++j)
            jmt_param_value(section, "java.lang.Boolean", bools[j], vals[j] ? "true" : "false");
        jmt_param_value(section, "java.lang.Integer", "numClasses",
                        jmt_int(static_cast<double>(sn_.nclasses)));
    }

    // -- model-level blocks -------------------------------------------------

    /**
     * Port of `saveMetrics`: one `<measure>` per enabled (station, class, kind),
     * plus the FCR and cache-hit-rate measures.
     *
     * `Residence Time` is NOT requested, as in the reference: JMT's definition
     * of it disagrees with LINE's on class-switching models, and the residence
     * time is recomputed from the response time and the visits instead.
     */
    void save_metrics(xml::Element& sim) {
        const JmtMetricKind kinds[5] = {JmtMetricKind::QLen, JmtMetricKind::Util,
                                        JmtMetricKind::RespT, JmtMetricKind::Tput,
                                        JmtMetricKind::ArvR};
        for (int k = 0; k < 5; ++k) save_metric(sim, kinds[k]);
        // Tardiness needs a class deadline, which this port's struct does not
        // carry; the measures are therefore not requested. See save_classes.
        save_fcr_metrics(sim);
        save_cache_hit_rate_metrics(sim);
    }

    void save_metric(xml::Element& sim, JmtMetricKind kind) {
        for (std::size_t ist = 1; ist <= sn_.nstations; ++ist)
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                if (!keep_[r - 1]) continue;
                if (!jmt_metric_enabled(sn_, kind, ist, r, cacheclass_)) continue;
                xml::Element& m = sim.add_child("measure");
                m.set_attr("alpha", jmt_sig2(1.0 - opt_.sim_conf_int));
                m.set_attr("name", std::string("Performance_") +
                                       jmt_int(static_cast<double>(ist)));
                m.set_attr("nodeType", "station");
                m.set_attr("precision", jmt_sig2(opt_.sim_max_rel_err));
                m.set_attr("referenceNode", nname(sn_.station_to_node[ist - 1]));
                m.set_attr("referenceUserClass", cname(r));
                m.set_attr("type", jmt_metric_text(kind));
                m.set_attr("verbose", "false");
            }
    }

    /** Port of `saveCacheHitRateMetrics`. */
    void save_cache_hit_rate_metrics(xml::Element& sim) {
        for (const auto& kv : sn_.nodeparam) {
            const std::size_t ind = kv.first;
            if (sn_.nodes[ind - 1].nodetype != NodeType::Cache) continue;
            const qn::CacheParam<T>& cp = kv.second;
            for (std::size_t r = 0; r < cp.hitclass.size(); ++r) {
                if (cp.hitclass[r] == 0) continue;
                xml::Element& m = sim.add_child("measure");
                m.set_attr("alpha", jmt_sig2(1.0 - opt_.sim_conf_int));
                m.set_attr("name", "CacheHitRate_" + nname(ind) + "_" + cname(cp.hitclass[r]));
                m.set_attr("nodeType", "station");
                m.set_attr("precision", jmt_sig2(opt_.sim_max_rel_err));
                m.set_attr("referenceNode", nname(ind));
                m.set_attr("referenceUserClass", cname(cp.hitclass[r]));
                m.set_attr("type", "Cache Hit Rate");
                m.set_attr("verbose", "false");
            }
        }
    }

    /** Port of `saveFCRMetrics`: six aggregate measures per region. */
    void save_fcr_metrics(xml::Element& sim) {
        static const char* kinds[6] = {"Number of Customers", "Response Time", "Residence Time",
                                       "Throughput",          "FCR Capacity",  "FCR Memory"};
        int counter = 0;
        for (std::size_t f = 1; f <= sn_.regions.size(); ++f) {
            const std::string fcr = "FCRegion" + jmt_int(static_cast<double>(f));
            for (int k = 0; k < 6; ++k) {
                std::string flat(kinds[k]);
                flat.erase(std::remove(flat.begin(), flat.end(), ' '), flat.end());
                xml::Element& m = sim.add_child("measure");
                m.set_attr("alpha", jmt_sig2(1.0 - opt_.sim_conf_int));
                m.set_attr("name", "FCR_" + fcr + "_" + flat + "_" +
                                       jmt_int(static_cast<double>(counter)));
                m.set_attr("nodeType", "region");
                m.set_attr("precision", jmt_sig2(opt_.sim_max_rel_err));
                m.set_attr("referenceNode", fcr);
                m.set_attr("referenceUserClass", "");
                m.set_attr("type", kinds[k]);
                m.set_attr("verbose", "false");
                ++counter;
            }
        }
    }

    /** Port of `saveLinks`: one `<connection>` per linked ordered pair. */
    void save_links(xml::Element& sim) {
        // Column-major, as MATLAB's `find` on the connection matrix returns.
        for (std::size_t j = 1; j <= sn_.nodes.size(); ++j)
            for (std::size_t i = 1; i <= sn_.nodes.size(); ++i) {
                if (!conn_[i - 1][j - 1]) continue;
                xml::Element& c = sim.add_child("connection");
                c.set_attr("source", nname(i));
                c.set_attr("target", nname(j));
            }
    }

    /**
     * Port of the `preload` block of `writeJSIM`: where the jobs start.
     *
     * A Source has an infinite reservoir and a Join holds no jobs, so neither
     * takes a preload. Every other station reports its per-class count, and a
     * closed class with no customers anywhere is omitted -- JMT reads an
     * omitted class as zero, and writing the zero would also declare the class
     * present at a station it never visits.
     */
    void save_preload(xml::Element& sim) {
        // Collected before anything is emitted: the block is omitted entirely
        // when no station takes a preload, and an element already appended to
        // the document cannot be withdrawn.
        struct Row {
            std::size_t node;
            std::vector<std::pair<std::size_t, double>> pops;  // (class, count)
        };
        std::vector<Row> rows;
        for (std::size_t ist = 1; ist <= sn_.nstations; ++ist) {
            const std::size_t ind = sn_.station_to_node[ist - 1];
            const NodeType ty = sn_.nodes[ind - 1].nodetype;
            if (ty == NodeType::Source || ty == NodeType::Join) continue;
            const std::vector<double> nir = initial_marginal(ist);
            Row row;
            row.node = ind;
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                const double n = sn_.classes[r - 1].population;
                if (std::isfinite(n) && n == 0.0 && nir[r - 1] == 0.0) continue;
                row.pops.emplace_back(r, nir[r - 1]);
            }
            if (!row.pops.empty()) rows.push_back(row);
        }
        if (rows.empty()) return;
        xml::Element& preload = sim.add_child("preload");
        for (const Row& row : rows) {
            xml::Element& st = preload.add_child("stationPopulations");
            st.set_attr("stationName", nname(row.node));
            for (const auto& pc : row.pops) {
                xml::Element& cp = st.add_child("classPopulation");
                cp.set_attr("population", jmt_int(pc.second));
                cp.set_attr("refClass", cname(pc.first));
            }
        }
    }

    /**
     * The per-class job count a station starts with.
     *
     * A Place's declared initial marking is its token count; every other
     * station takes `model.initDefault`, which puts each closed class's whole
     * population at its reference station and nothing anywhere else. Open
     * classes start empty.
     */
    std::vector<double> initial_marginal(std::size_t ist) const {
        std::vector<double> nir(sn_.nclasses, 0.0);
        const std::size_t ind = sn_.station_to_node[ist - 1];
        const auto im = sn_.initmarking.find(ind);
        if (im != sn_.initmarking.end()) {
            for (std::size_t r = 0; r < sn_.nclasses && r < im->second.size(); ++r)
                nir[r] = d(im->second[r]);
            return nir;
        }
        for (std::size_t r = 0; r < sn_.nclasses; ++r) {
            const double n = sn_.classes[r].population;
            if (std::isfinite(n) && sn_.classes[r].refstat == ist) nir[r] = n;
        }
        return nir;
    }

    // -- Place (Storage / Linkage) ------------------------------------------

    /** Port of `saveTotalCapacity`: the token bound of a place, -1 unbounded. */
    void save_total_capacity(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        const std::string v =
            (ist == 0 || !std::isfinite(sn_.cap[ist - 1])) ? "-1" : jmt_int(sn_.cap[ist - 1]);
        jmt_param_value(section, "java.lang.Integer", "totalCapacity", v);
    }

    /** Port of `savePlaceCapacities`: the per-colour token bound. */
    void save_place_capacities(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "java.lang.Integer", "capacities", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            const double c = ist == 0 ? std::numeric_limits<double>::infinity()
                                      : sn_.classcap[ist - 1][r - 1];
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Integer");
            sp.set_attr("name", "capacity");
            sp.add_text_child("value", std::isfinite(c) ? jmt_int(c) : "-1");
        }
    }

    /**
     * Port of `saveDropRule`. It differs from `saveDropStrategy` only in the
     * parameter names JMT's Storage section expects (`dropRules`/`dropRule`).
     */
    void save_drop_rule(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "java.lang.String", "dropRules", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.String");
            sp.set_attr("name", "dropRule");
            sp.add_text_child("value", drop_strategy_text(ist, r));
        }
    }

    /**
     * Port of `savePutStrategies` (plural), the Storage section's insertion
     * rule. It offers only the three orders a place can hold tokens in, which
     * is why it is a separate handler from `savePutStrategy`.
     */
    void save_put_strategies(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(section, "jmt.engine.NetStrategies.QueuePutStrategy",
                                    "QueuePutStrategy", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            const char* nm = "TailStrategy";
            if (ist != 0) {
                if (sn_.stations[ist - 1].sched == SchedStrategy::SIRO) nm = "RandStrategy";
                else if (sn_.stations[ist - 1].sched == SchedStrategy::LCFS) nm = "HeadStrategy";
            }
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath",
                        std::string("jmt.engine.NetStrategies.QueuePutStrategies.") + nm);
            sp.set_attr("name", nm);
        }
    }

    // -- Transition (Enabling / Timing / Firing) ----------------------------

    /** The transition parameters of node `ind`; refuses a node that has none. */
    const qn::TransitionParam<T>& transparam(std::size_t ind) const {
        const auto it = sn_.transparam.find(ind);
        if (it == sn_.transparam.end())
            throw InputError("SolverJMT: transition '" + nname(ind) + "' carries no mode table");
        return it->second;
    }

    /**
     * Port of `saveEnablingConditions`.
     *
     * A place appears in the vector only when it carries a POSITIVE enabling or
     * inhibiting entry for the mode IN SOME CLASS. An enabling entry of infinity
     * is JMT's -1, "any number of tokens".
     *
     * THE ENTRIES ARE PER CLASS, which is what JSIM's own vector format is: one
     * `enablingEntry` per class, in class order. Writing the same number under
     * every class -- what this writer did while the struct had no class
     * dimension -- exported a net that demands every colour on every arc, and
     * that is why a multiclass net had to be refused here.
     */
    void save_enabling_conditions(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        xml::Element& p =
            jmt_param(section, "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix",
                      "enablingConditions", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            xml::Element& cond = p.add_child("subParameter");
            cond.set_attr("classPath",
                          "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
            cond.set_attr("name", "enablingCondition");
            xml::Element& vecs = jmt_param_sub(
                cond, "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector",
                "enablingVectors");
            for (std::size_t k = 1; k <= sn_.nodes.size(); ++k) {
                bool relevant = false;
                for (std::size_t r = 1; r <= sn_.nclasses && !relevant; ++r) {
                    const double en = arc(tp.enabling, m, k, r);
                    const double in = arc(tp.inhibiting, m, k, r);
                    relevant = (std::isfinite(en) && en > 0.0) ||
                               (std::isfinite(in) && in > 0.0);
                }
                if (!relevant) continue;
                xml::Element& vec = vecs.add_child("subParameter");
                vec.set_attr("classPath",
                             "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
                vec.set_attr("name", "enablingVector");
                jmt_scalar(vec, "java.lang.String", "stationName", nname(k));
                xml::Element& entries = jmt_param_sub(vec, "java.lang.Integer", "enablingEntries");
                for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                    const double en = arc(tp.enabling, m, k, r);
                    jmt_ref_class(entries, cname(r));
                    xml::Element& e = entries.add_child("subParameter");
                    e.set_attr("classPath", "java.lang.Integer");
                    e.set_attr("name", "enablingEntry");
                    e.add_text_child("value", std::isfinite(en) ? jmt_int(en) : "-1");
                }
            }
        }
    }

    /**
     * Port of `saveInhibitingConditions`.
     *
     * The vectors cover the transition's INPUT places only, and an infinite
     * entry -- "no inhibitor arc" -- is written as 0, which is what JMT reads as
     * "never inhibits". Note that the reference's sentinel is inverted between
     * the two handlers: infinity is -1 in the enabling block and 0 here.
     */
    void save_inhibiting_conditions(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        const std::vector<std::size_t> inputs = inputs_of(ind);
        xml::Element& p =
            jmt_param(section, "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix",
                      "inhibitingConditions", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            xml::Element& cond = p.add_child("subParameter");
            cond.set_attr("classPath",
                          "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
            cond.set_attr("name", "inhibitingCondition");
            xml::Element& vecs = jmt_param_sub(
                cond, "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector",
                "inhibitingVectors");
            for (std::size_t k : inputs) {
                xml::Element& vec = vecs.add_child("subParameter");
                vec.set_attr("classPath",
                             "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
                vec.set_attr("name", "inhibitingVector");
                jmt_scalar(vec, "java.lang.String", "stationName", nname(k));
                xml::Element& entries = jmt_param_sub(vec, "java.lang.Integer", "inhibitingEntries");
                for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                    const double in = arc(tp.inhibiting, m, k, r);
                    jmt_ref_class(entries, cname(r));
                    xml::Element& e = entries.add_child("subParameter");
                    e.set_attr("classPath", "java.lang.Integer");
                    e.set_attr("name", "inhibitingEntry");
                    e.add_text_child("value", std::isfinite(in) ? jmt_int(in) : "0");
                }
            }
        }
    }

    /** Port of `saveFiringOutcomes`: the tokens a firing puts in each output. */
    void save_firing_outcomes(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        const std::vector<std::size_t> outs = outputs_of(ind);
        xml::Element& p =
            jmt_param(section, "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix",
                      "firingOutcomes", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            xml::Element& out = p.add_child("subParameter");
            out.set_attr("classPath",
                         "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
            out.set_attr("name", "firingOutcome");
            xml::Element& vecs = jmt_param_sub(
                out, "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector",
                "firingVectors");
            for (std::size_t k : outs) {
                xml::Element& vec = vecs.add_child("subParameter");
                vec.set_attr("classPath",
                             "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
                vec.set_attr("name", "firingVector");
                jmt_scalar(vec, "java.lang.String", "stationName", nname(k));
                xml::Element& entries = jmt_param_sub(vec, "java.lang.Integer", "firingEntries");
                for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                    const double f = arc(tp.firing, m, k, r);
                    jmt_ref_class(entries, cname(r));
                    xml::Element& e = entries.add_child("subParameter");
                    e.set_attr("classPath", "java.lang.Integer");
                    e.set_attr("name", "firingEntry");
                    e.add_text_child("value", jmt_int(f));
                }
            }
        }
    }

    /** `enabling`/`inhibiting`/`firing` of mode m at 1-based node k, class r. */
    double arc(const std::vector<Matrix<T>>& tab, std::size_t m, std::size_t k,
               std::size_t r) const {
        if (m >= tab.size() || k == 0 || k > tab[m].rows() || r == 0 || r > tab[m].cols())
            return 0.0;
        return d(tab[m](k - 1, r - 1));
    }

    /** Port of `saveModeNames`. */
    void save_mode_names(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        xml::Element& p = jmt_param(section, "java.lang.String", "modeNames", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.String");
            sp.set_attr("name", "modeName");
            sp.add_text_child("value", m < tp.modenames.size() ? tp.modenames[m] : std::string());
        }
    }

    /** Port of `saveNumbersOfServers`: the firing concurrency of each mode. */
    void save_numbers_of_servers(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        xml::Element& p = jmt_param(section, "java.lang.Integer", "numbersOfServers", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            const double ns = m < tp.nmodeservers.size() ? tp.nmodeservers[m] : 1.0;
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Integer");
            sp.set_attr("name", "numberOfServers");
            sp.add_text_child("value", std::isfinite(ns) ? jmt_int(ns) : "-1");
        }
    }

    /**
     * Port of `saveTimingStrategies`.
     *
     * An immediate mode is a ZeroServiceTimeStrategy and carries no
     * distribution; a timed one carries its firing process through the shared
     * emitter, under JMT's `timingStrategy` name rather than the
     * `ServiceTimeStrategy` name a queue's service uses.
     */
    void save_timing_strategies(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        xml::Element& p = jmt_param(section, "jmt.engine.NetStrategies.ServiceStrategy",
                                    "timingStrategies", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            xml::Element& sp = p.add_child("subParameter");
            const bool immediate = m < tp.timing.size() &&
                                   tp.timing[m] == lang::TimingStrategy::IMMEDIATE;
            if (immediate) {
                sp.set_attr("classPath",
                            "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                sp.set_attr("name", "ZeroServiceTimeStrategy");
                continue;
            }
            sp.set_attr("classPath",
                        "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            sp.set_attr("name", "timingStrategy");
            if (m >= tp.firingproc.size())
                throw InputError("SolverJMT: timed mode '" +
                                 (m < tp.modenames.size() ? tp.modenames[m] : std::string()) +
                                 "' of transition '" + nname(ind) + "' has no firing process");
            jmt_append_distribution(sp, jmt_dist_view(tp.firingproc[m]),
                                    "SolverJMT (transition firing)");
        }
    }

    /** Port of `saveFiringPriorities`. */
    void save_firing_priorities(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        xml::Element& p = jmt_param(section, "java.lang.Integer", "firingPriorities", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            const double v = m < tp.firingprio.size() ? tp.firingprio[m] : 1.0;
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Integer");
            sp.set_attr("name", "firingPriority");
            sp.add_text_child("value", std::isfinite(v) ? jmt_int(v) : "-1");
        }
    }

    /**
     * Port of `saveFiringWeights`.
     *
     * THE REFERENCE PRINTS THE WEIGHT AS AN INTEGER (`int2str`) although the
     * parameter is a `java.lang.Double`, so a weight of 0.3 exports as 0 and
     * the mode never wins a race it should sometimes win. That is a defect and
     * not a convention -- no other weight in the file is rounded -- so this
     * port writes the weight itself.
     */
    void save_firing_weights(xml::Element& section, std::size_t ind) {
        const qn::TransitionParam<T>& tp = transparam(ind);
        xml::Element& p = jmt_param(section, "java.lang.Double", "firingWeights", true);
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            const double v = m < tp.fireweight.size() ? d(tp.fireweight[m]) : 1.0;
            xml::Element& sp = p.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Double");
            sp.set_attr("name", "firingWeight");
            sp.add_text_child("value", std::isfinite(v) ? jmt_fmt(v) : "-1");
        }
    }

    // -- blocking regions ---------------------------------------------------

    /**
     * Port of `jmtClassCapCon`: the per-(station, class) capacities that are a
     * REAL constraint and must therefore reach JMT as a blocking region.
     *
     * A capacity is real only when it binds: not at a Source or a Place, not on
     * a pair the class never visits, and not when it already exceeds the
     * station's own capacity or the population of the class's chain -- such a
     * bound can never be reached, and exporting it as a region would add a
     * region JMT then reports measures for.
     */
    std::vector<std::vector<double>> class_cap_constraints() const {
        const double inf = std::numeric_limits<double>::infinity();
        std::vector<std::vector<double>> con(sn_.nstations, std::vector<double>(sn_.nclasses, inf));
        if (sn_.classcap.empty()) return con;
        std::vector<double> chainpop(sn_.nclasses, inf);
        for (std::size_t c = 0; c < sn_.inchain.size(); ++c) {
            double tot = 0.0;
            bool open = false;
            for (std::size_t r : sn_.inchain[c]) {
                const double n = sn_.classes[r - 1].population;
                if (std::isfinite(n)) tot += n; else open = true;
            }
            for (std::size_t r : sn_.inchain[c]) chainpop[r - 1] = open ? inf : tot;
        }
        for (std::size_t ist = 1; ist <= sn_.nstations; ++ist) {
            const NodeType ty = sn_.nodes[sn_.station_to_node[ist - 1] - 1].nodetype;
            if (ty == NodeType::Source || ty == NodeType::Place) continue;
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                if (sn_.disabled[ist - 1][r - 1]) continue;
                const double cc = sn_.classcap[ist - 1][r - 1];
                if (std::isfinite(cc) && cc < 2147483647.0 &&
                    cc < std::min(sn_.cap[ist - 1], chainpop[r - 1]))
                    con[ist - 1][r - 1] = cc;
            }
        }
        return con;
    }

    /**
     * Port of `jmtClassCapAssert`: the two per-class capacities JMT cannot
     * reproduce, refused by name rather than exported as something else.
     *
     * A CLOSED class is the first. LINE holds a blocked closed job at its
     * upstream station; JMT's blocking region parks it in the region's input
     * station instead, which frees the upstream server and loses the job from
     * the population count -- a different model, not a different rounding.
     * A non-loss drop strategy is the second: a region can only drop or defer.
     */
    void assert_class_cap_exportable(std::size_t ist, std::size_t r) const {
        const std::string sname = nname(sn_.station_to_node[ist - 1]);
        if (std::isfinite(sn_.classes[r - 1].population))
            throw UnsupportedError(
                "SolverJMT: station '" + sname + "' carries a finite capacity " +
                jmt_int(sn_.classcap[ist - 1][r - 1]) + " for the closed class '" + cname(r) +
                "'. LINE holds a blocked closed job at its upstream station, whereas JMT can "
                "only express a per-class capacity as a blocking region, which parks the job in "
                "the region input station instead, freeing the upstream server and losing it "
                "from the population count. Use SolverCTMC, SolverSSA or SolverLDES, or express "
                "the limit as the station capacity");
        const DropStrategy dr = sn_.droprule[ist - 1][r - 1];
        if (dr == DropStrategy::BAS || dr == DropStrategy::BBS || dr == DropStrategy::RSRD ||
            dr == DropStrategy::RETRIAL || dr == DropStrategy::RETRIAL_WITH_LIMIT)
            throw UnsupportedError(
                "SolverJMT: station '" + sname + "' applies drop strategy '" + jmt_drop_text(dr) +
                "' to class '" + cname(r) +
                "' and also carries a finite capacity for it. JMT exports a per-class capacity "
                "as a blocking region, which can only drop or defer an arrival. Remove the "
                "per-class capacity, or use the station capacity, which is exported with its "
                "drop strategy");
    }

    /**
     * Port of `saveRegions`: three kinds of `<blockingRegion>`, in this order.
     *
     * 1. ONE PER LPS STATION. Limited processor sharing is a cap on the number
     *    IN SERVICE, which JMT has no station-level field for; the single-node
     *    region is how the limit is expressed, and it is why
     *    `save_number_of_servers` exports an LPS station as a single server.
     * 2. The model's declared finite capacity regions.
     * 3. ONE PER STATION whose per-class capacity is a real constraint and that
     *    no region above already covers, since JMT allows a node to belong to
     *    at most one region.
     */
    void save_regions(xml::Element& sim) {
        const std::vector<std::vector<double>> con = class_cap_constraints();
        std::vector<bool> covered(sn_.nstations, false);
        std::size_t lps_idx = sn_.regions.size();

        for (std::size_t ist = 1; ist <= sn_.nstations; ++ist) {
            if (sn_.stations[ist - 1].sched != SchedStrategy::LPS) continue;
            ++lps_idx;
            const std::size_t ind = sn_.station_to_node[ist - 1];
            const double limit = sn_.stations[ist - 1].schedparam.empty()
                                     ? 1.0
                                     : d(sn_.stations[ist - 1].schedparam[0]);
            xml::Element& br = sim.add_child("blockingRegion");
            br.set_attr("name", "LPSRegion" + jmt_int(static_cast<double>(lps_idx)));
            br.set_attr("type", "default");
            br.add_child("regionNode").set_attr("nodeName", nname(ind));
            br.add_child("globalConstraint").set_attr("maxJobs", jmt_num(limit));
            br.add_child("globalMemoryConstraint").set_attr("maxMemory", "-1");
            covered[ist - 1] = true;
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                if (!std::isfinite(con[ist - 1][r - 1])) continue;
                assert_class_cap_exportable(ist, r);  // rejects the non-loss cases first
                throw UnsupportedError(
                    "SolverJMT: station '" + nname(ind) +
                    "' has both LPS scheduling and a finite capacity for the open class '" +
                    cname(r) +
                    "'. JMT expresses both through a single blocking region, which admits only "
                    "one drop rule per class, but LPS requires blocking while the open-class "
                    "capacity requires dropping. Remove the per-class capacity or use a non-LPS "
                    "scheduling strategy");
            }
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                xml::Element& dr = br.add_child("dropRules");
                dr.set_attr("jobClass", cname(r));
                dr.set_attr("dropThisClass", "false");
            }
        }

        for (std::size_t f = 1; f <= sn_.regions.size(); ++f) {
            const typename qn::NetworkStruct<T>::Region& rg = sn_.regions[f - 1];
            std::vector<std::size_t> members;
            for (std::size_t ist = 1; ist <= sn_.nstations; ++ist)
                if (ist <= rg.members.size() && rg.members[ist - 1]) {
                    members.push_back(ist);
                    covered[ist - 1] = true;
                }
            std::vector<double> region_class_cap(sn_.nclasses,
                                                 std::numeric_limits<double>::infinity());
            for (std::size_t ist : members)
                for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                    if (!std::isfinite(con[ist - 1][r - 1])) continue;
                    if (members.size() > 1)
                        throw UnsupportedError(
                            "SolverJMT: station '" + nname(sn_.station_to_node[ist - 1]) +
                            "' carries a finite capacity for class '" + cname(r) +
                            "' and also belongs to the multi-station region FCRegion" +
                            jmt_int(static_cast<double>(f)) +
                            ". JMT constrains a blocking region as a whole and allows a node to "
                            "belong to only one region, so a per-station class capacity cannot "
                            "be expressed alongside it");
                    assert_class_cap_exportable(ist, r);
                    if (!(r <= rg.rule.size() && rg.rule[r - 1] == DropStrategy::DROP))
                        throw UnsupportedError(
                            "SolverJMT: station '" + nname(sn_.station_to_node[ist - 1]) +
                            "' carries a finite capacity for class '" + cname(r) +
                            "' and also belongs to region FCRegion" +
                            jmt_int(static_cast<double>(f)) +
                            ", whose drop rule for that class is not DROP. A JMT blocking "
                            "region admits a single drop rule per class, shared by all of its "
                            "constraints, and the per-class capacity of an open class is a loss "
                            "constraint");
                    region_class_cap[r - 1] = con[ist - 1][r - 1];
                }

            xml::Element& br = sim.add_child("blockingRegion");
            br.set_attr("name", "FCRegion" + jmt_int(static_cast<double>(f)));
            br.set_attr("type", "default");
            for (std::size_t ist : members)
                br.add_child("regionNode")
                    .set_attr("nodeName", nname(sn_.station_to_node[ist - 1]));
            // The global caps are stored per MEMBER station and are the same at
            // each; the first member therefore carries the region's own values.
            const double gmax =
                members.empty() ? -1.0 : rg.cap[members[0] - 1][sn_.nclasses];
            const double gmem = members.empty() ? -1.0 : rg.maxmem[members[0] - 1];
            br.add_child("globalConstraint").set_attr("maxJobs", jmt_num(gmax));
            br.add_child("globalMemoryConstraint").set_attr("maxMemory", jmt_num(gmem));
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                double cmax = members.empty() ? -1.0 : rg.cap[members[0] - 1][r - 1];
                // -1 is the struct's "unbounded"; the per-station capacity then
                // stands alone, and otherwise the tighter of the two binds.
                if (cmax == -1.0)
                    cmax = std::isfinite(region_class_cap[r - 1]) ? region_class_cap[r - 1] : -1.0;
                else if (std::isfinite(region_class_cap[r - 1]))
                    cmax = std::min(cmax, region_class_cap[r - 1]);
                if (cmax == -1.0 || !std::isfinite(cmax)) continue;
                xml::Element& cc = br.add_child("classConstraint");
                cc.set_attr("jobClass", cname(r));
                cc.set_attr("maxJobsPerClass", jmt_num(cmax));
            }
            // NO `classMemoryConstraint` IS EMITTED. `add_region` already folds
            // the per-class memory budget into the per-class job cap by
            // dividing it by the class size, so the constraint above carries it;
            // emitting it again would apply the same budget twice.
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                xml::Element& dr = br.add_child("dropRules");
                dr.set_attr("jobClass", cname(r));
                dr.set_attr("dropThisClass",
                            (r <= rg.rule.size() && rg.rule[r - 1] == DropStrategy::DROP)
                                ? "true"
                                : "false");
            }
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                if (r > rg.weight.size()) continue;
                const double w = d(rg.weight[r - 1]);
                if (w == 1.0) continue;
                xml::Element& cw = br.add_child("classWeight");
                cw.set_attr("jobClass", cname(r));
                cw.set_attr("weight", jmt_num(w));
            }
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                if (r > rg.size.size()) continue;
                const double s = d(rg.size[r - 1]);
                if (s == 1.0) continue;
                xml::Element& cs = br.add_child("classSize");
                cs.set_attr("jobClass", cname(r));
                cs.set_attr("size", jmt_num(s));
            }
        }

        std::size_t cap_idx = 0;
        for (std::size_t ist = 1; ist <= sn_.nstations; ++ist) {
            if (covered[ist - 1]) continue;
            bool any = false;
            for (std::size_t r = 1; r <= sn_.nclasses; ++r)
                if (std::isfinite(con[ist - 1][r - 1])) any = true;
            if (!any) continue;
            ++cap_idx;
            for (std::size_t r = 1; r <= sn_.nclasses; ++r)
                if (std::isfinite(con[ist - 1][r - 1])) assert_class_cap_exportable(ist, r);
            xml::Element& br = sim.add_child("blockingRegion");
            br.set_attr("name", "ClassCapRegion" + jmt_int(static_cast<double>(cap_idx)));
            br.set_attr("type", "default");
            br.add_child("regionNode")
                .set_attr("nodeName", nname(sn_.station_to_node[ist - 1]));
            br.add_child("globalConstraint").set_attr("maxJobs", "-1");
            br.add_child("globalMemoryConstraint").set_attr("maxMemory", "-1");
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                if (!std::isfinite(con[ist - 1][r - 1])) continue;
                xml::Element& cc = br.add_child("classConstraint");
                cc.set_attr("jobClass", cname(r));
                cc.set_attr("maxJobsPerClass", jmt_num(con[ist - 1][r - 1]));
            }
            for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
                if (!std::isfinite(con[ist - 1][r - 1])) continue;
                xml::Element& dr = br.add_child("dropRules");
                dr.set_attr("jobClass", cname(r));
                dr.set_attr("dropThisClass", "true");
            }
        }
    }

    // -- Cache --------------------------------------------------------------

    /**
     * Port of `saveCacheStrategy`.
     *
     * THE LIST-TO-LIST MATRIX IS DERIVED FROM THE POLICY, not from `accost`:
     * LRU promotes an item one list up on a hit and the last list keeps it,
     * every other policy keeps it where it is. That is what the reference
     * emits, and JMT's own replacement object does the rest.
     *
     * SFIFO is exported as JMT's FIFO cache and every policy JMT has no object
     * for -- HLRU, CLIMB, QLRU -- is refused by name rather than falling back
     * to LRU as the reference's `otherwise` does: a q-LRU cache silently
     * simulated as LRU returns a hit rate that is wrong by the admission
     * probability, with nothing in the output to show for it.
     */
    void save_cache_strategy(xml::Element& section, std::size_t ind) {
        const auto it = sn_.nodeparam.find(ind);
        if (it == sn_.nodeparam.end())
            throw InputError("SolverJMT: cache node '" + nname(ind) + "' carries no parameters");
        const qn::CacheParam<T>& cp = it->second;
        const std::size_t K = sn_.nclasses;

        jmt_param_value(section, "java.lang.Integer", "maxItems",
                        jmt_int(static_cast<double>(cp.nitems)));
        xml::Element& cap = jmt_param(section, "java.lang.Integer", "cacheCapacity", true);
        for (std::size_t l = 0; l < cp.itemcap.size(); ++l) {
            xml::Element& sp = cap.add_child("subParameter");
            sp.set_attr("classPath", "java.lang.Integer");
            sp.set_attr("name", "capacity");
            sp.add_text_child("value", jmt_int(static_cast<double>(cp.itemcap[l])));
        }

        const std::size_t nlev = cp.itemcap.size();
        const bool lru = cp.replacestrat == lang::ReplacementStrategy::LRU;
        xml::Element& mat = jmt_param(section, "java.lang.Object", "matrix", true);
        for (std::size_t a = 1; a <= nlev; ++a) {
            xml::Element& row = mat.add_child("subParameter");
            row.set_attr("array", "true");
            row.set_attr("classPath", "java.lang.Float");
            row.set_attr("name", "row");
            for (std::size_t b = 1; b <= nlev; ++b) {
                const bool one = lru ? (a < nlev ? b == a + 1 : b == nlev) : a == b;
                xml::Element& cell = row.add_child("subParameter");
                cell.set_attr("classPath", "java.lang.Float");
                cell.set_attr("name", "cell");
                cell.add_text_child("value", one ? "1.0" : "0.0");
            }
        }

        xml::Element& jc = jmt_param(section, "jmt.engine.QueueNet.JobClass", "jobClasses", true);
        for (std::size_t r = 1; r <= K; ++r) {
            const bool used = (r <= cp.hitclass.size() && cp.hitclass[r - 1] > 0) ||
                              (r <= cp.missclass.size() && cp.missclass[r - 1] > 0);
            if (!used) continue;
            xml::Element& sp = jc.add_child("subParameter");
            sp.set_attr("classPath", "jmt.engine.QueueNet.JobClass");
            sp.set_attr("name", "jobClass");
            sp.add_text_child("value", cname(r));
        }
        const char* switch_names[2] = {"hitClasses", "missClasses"};
        const char* entry_names[2] = {"hitClass", "missClass"};
        for (int w = 0; w < 2; ++w) {
            const std::vector<std::size_t>& tab = w == 0 ? cp.hitclass : cp.missclass;
            xml::Element& p =
                jmt_param(section, "jmt.engine.QueueNet.JobClass", switch_names[w], true);
            for (std::size_t r = 1; r <= K && r <= tab.size(); ++r) {
                if (tab[r - 1] == 0) continue;
                xml::Element& sp = p.add_child("subParameter");
                sp.set_attr("classPath", "jmt.engine.QueueNet.JobClass");
                sp.set_attr("name", entry_names[w]);
                sp.add_text_child("value", cname(tab[r - 1]));
            }
        }

        const char* policy = nullptr;
        switch (cp.replacestrat) {
            case lang::ReplacementStrategy::LRU:
                policy = "jmt.engine.NetStrategies.CacheStrategies.LRUCache";
                break;
            case lang::ReplacementStrategy::FIFO:
            case lang::ReplacementStrategy::SFIFO:
                policy = "jmt.engine.NetStrategies.CacheStrategies.FIFOCache";
                break;
            case lang::ReplacementStrategy::RR:
                policy = "jmt.engine.NetStrategies.CacheStrategies.RandomCache";
                break;
            default:
                throw UnsupportedError(
                    "SolverJMT: cache '" + nname(ind) +
                    "' uses a replacement policy JMT has no cache object for (HLRU, CLIMB and "
                    "QLRU); use SolverLDES, which simulates them directly");
        }
        xml::Element& rp = section.add_child("parameter");
        rp.set_attr("classPath", policy);
        rp.set_attr("name", "replacePolicy");

        /*
         * The popularity is PARAMETRIC in JMT: a Zipf exponent or a uniform
         * range, never a pmf. `CacheParam::preadkind` is what the reader kept
         * for exactly this, and a class whose pmf was supplied directly exports
         * as `null` -- the reference's own behaviour for anything that is not a
         * Zipf or a DiscreteSampler.
         */
        xml::Element& pop = jmt_param(
            section, "jmt.engine.random.discrete.DiscreteDistribution", "popularity", true);
        for (std::size_t r = 1; r <= K; ++r) {
            jmt_ref_class(pop, cname(r));
            const bool reads = r <= cp.pread.size() && !cp.pread[r - 1].empty();
            const typename qn::CacheParam<T>::Popularity kind =
                r <= cp.preadkind.size() ? cp.preadkind[r - 1]
                                         : typename qn::CacheParam<T>::Popularity();
            xml::Element& sp = pop.add_child("subParameter");
            if (reads && kind.type == lang::ProcessType::ZIPF) {
                sp.set_attr("classPath", "jmt.engine.random.discrete.Zipf");
                sp.set_attr("name", "popularity");
                jmt_scalar(sp, "java.lang.Double", "alpha", jmt_fmt(kind.s));
                jmt_scalar(sp, "java.lang.Integer", "numberOfElements",
                           jmt_int(static_cast<double>(kind.n)));
            } else if (reads && kind.type == lang::ProcessType::DISCRETESAMPLER) {
                sp.set_attr("classPath", "jmt.engine.random.discrete.Uniform");
                sp.set_attr("name", "popularity");
                jmt_scalar(sp, "java.lang.Integer", "min", "1");
                jmt_scalar(sp, "java.lang.Integer", "max",
                           jmt_int(static_cast<double>(kind.n == 0 ? cp.nitems : kind.n)));
            } else {
                sp.set_attr("classPath", "jmt.engine.random.discrete.DiscreteDistribution");
                sp.set_attr("name", "null");
                sp.add_text_child("value", "null");
            }
        }
    }

    // -- impatience, retrials and switchover --------------------------------

    /**
     * Port of `saveImpatience`: the per-class Balking or Reneging strategy.
     *
     * BALKING AND RENEGING ARE MUTUALLY EXCLUSIVE per class here, as in the
     * reference: JMT's Impatience array holds one strategy per class, and
     * balking is checked first. A class with neither gets a Reneging strategy
     * whose body is the literal `null`.
     */
    void save_impatience(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        xml::Element& p = jmt_param(
            section, "jmt.engine.NetStrategies.ImpatienceStrategies.Impatience", "Impatience", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& sp = p.add_child("subParameter");
            const qn::Station<T>* st = ist == 0 ? nullptr : &sn_.stations[ist - 1];
            const bool balks = st != nullptr && r <= st->balking.size() &&
                               st->balking[r - 1].strategy ==
                                   lang::BalkingStrategy::QUEUE_LENGTH;
            if (balks) {
                sp.set_attr("classPath", "jmt.engine.NetStrategies.ImpatienceStrategies.Balking");
                sp.set_attr("name", "Balking");
                save_balking_strategy(sp, st->balking[r - 1].thresholds, st->nservers);
                continue;
            }
            sp.set_attr("classPath", "jmt.engine.NetStrategies.ImpatienceStrategies.Reneging");
            sp.set_attr("name", "Reneging");
            const bool renegs = st != nullptr && r <= st->impatience.size() &&
                                st->impatience[r - 1] == lang::ImpatienceType::RENEGING &&
                                r <= st->patience.size() && !st->patience[r - 1].disabled;
            if (!renegs) {
                sp.add_text_child("value", "null");
                continue;
            }
            const JmtDistView<T> v = jmt_dist_view(st->patience[r - 1]);
            if (v.type == lang::ProcessType::UNIFORM) {
                // The reference exports a patience Uniform as [0, 2*mean],
                // matching the MEAN only, where the service-time emitter
                // matches the mean and the SCV. Kept: a patience law is
                // declared by its mean in every model that has one.
                xml::Element& dn = sp.add_child("subParameter");
                dn.set_attr("classPath", "jmt.engine.random.Uniform");
                dn.set_attr("name", "Uniform");
                xml::Element& par = sp.add_child("subParameter");
                par.set_attr("classPath", "jmt.engine.random.UniformPar");
                par.set_attr("name", "distrPar");
                jmt_scalar(par, "java.lang.Double", "min", "0.0");
                jmt_double(par, "max", v.rate > 0.0 ? 2.0 / v.rate : 0.0);
                continue;
            }
            jmt_append_distribution(sp, v, "SolverJMT (patience)");
        }
    }

    /**
     * Port of `saveBalkingStrategy`.
     *
     * THE THRESHOLDS ARE SHIFTED DOWN BY THE SERVER COUNT. LINE evaluates
     * balking against the TOTAL station population, JMT's Balking against the
     * number WAITING; in the regime where balking matters the servers are busy,
     * so the two differ by exactly S and JMT waiting w is LINE total w + S.
     * A gap between two intervals gets an explicit zero-probability breakpoint,
     * since JMT selects the LAST range with `from <= queueLength` and would
     * otherwise let a queue length outside every interval inherit its
     * neighbour's balking probability.
     */
    void save_balking_strategy(xml::Element& balking,
                               const std::vector<typename qn::Station<T>::BalkingThreshold>& th,
                               double nservers) {
        struct Bp {
            double from, prob;
        };
        std::vector<typename qn::Station<T>::BalkingThreshold> sorted = th;
        std::sort(sorted.begin(), sorted.end(),
                  [](const typename qn::Station<T>::BalkingThreshold& a,
                     const typename qn::Station<T>::BalkingThreshold& b) {
                      return a.min_jobs < b.min_jobs;
                  });
        const double S =
            (!std::isfinite(nservers) || nservers < 1.0) ? 1.0 : nservers;
        std::vector<Bp> bps;
        for (std::size_t i = 0; i < sorted.size(); ++i) {
            const double lo = sorted[i].min_jobs;
            // -1 is the wire's unbounded upper end; every other value is a
            // closed interval, which is why the gap test below is `hi + 1`.
            const double hi = sorted[i].max_jobs < 0.0
                                  ? std::numeric_limits<double>::infinity()
                                  : sorted[i].max_jobs;
            const double pr = d(sorted[i].probability);
            bps.push_back(Bp{std::max(0.0, lo - S), pr});
            if (!std::isfinite(hi)) continue;
            const double next_lo = i + 1 < sorted.size()
                                       ? sorted[i + 1].min_jobs
                                       : std::numeric_limits<double>::infinity();
            if (hi + 1.0 < next_lo) bps.push_back(Bp{std::max(0.0, hi + 1.0 - S), 0.0});
        }
        xml::Element& ld = balking.add_child("subParameter");
        ld.set_attr("classPath", "jmt.engine.NetStrategies.ServiceStrategies.LoadDependentStrategy");
        ld.set_attr("name", "LoadDependentStrategy");
        xml::Element& arr = jmt_param_sub(
            ld, "jmt.engine.NetStrategies.ServiceStrategies.LDParameter", "LDParameter");
        for (const Bp& b : bps) {
            xml::Element& rn = arr.add_child("subParameter");
            rn.set_attr("classPath", "jmt.engine.NetStrategies.ServiceStrategies.LDParameter");
            rn.set_attr("name", "LDParameter");
            jmt_scalar(rn, "java.lang.Integer", "from", jmt_int(b.from));
            // A DUMMY distribution: only the `function` string below is read
            // for balking, but the LDParameter shape requires the pair.
            xml::Element& dn = rn.add_child("subParameter");
            dn.set_attr("classPath", "jmt.engine.random.Exponential");
            dn.set_attr("name", "Exponential");
            xml::Element& par = rn.add_child("subParameter");
            par.set_attr("classPath", "jmt.engine.random.ExponentialPar");
            par.set_attr("name", "distrPar");
            jmt_scalar(par, "java.lang.Double", "lambda", "1.0");
            jmt_scalar(rn, "java.lang.String", "function", jmt_fmt(b.prob));
        }
        jmt_scalar(balking, "java.lang.Boolean", "priorityActivated", "false");
    }

    /**
     * Port of `saveRetrialDistributions`: the orbit delay of each class.
     *
     * A class with no orbit still gets an entry -- an Exp(1) placeholder -- so
     * that JMT's per-class array stays aligned; omitting it would shift every
     * later class's retrial delay onto the wrong class.
     */
    void save_retrial_distributions(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        const auto it = ist == 0 ? sn_.retrialparam.end() : sn_.retrialparam.find(ist);
        xml::Element& p = jmt_param(section, "jmt.engine.NetStrategies.ServiceStrategy",
                                    "retrialDistributions", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            if (!keep_[r - 1]) continue;
            jmt_ref_class(p, cname(r));
            xml::Element& sts = p.add_child("subParameter");
            sts.set_attr("classPath",
                         "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            sts.set_attr("name", "ServiceTimeStrategy");
            const bool has = it != sn_.retrialparam.end() &&
                             r <= it->second.retrial_proc.size() &&
                             !it->second.retrial_proc[r - 1].disabled;
            if (!has) {
                xml::Element& dn = sts.add_child("subParameter");
                dn.set_attr("classPath", "jmt.engine.random.Exponential");
                dn.set_attr("name", "Exponential");
                xml::Element& par = sts.add_child("subParameter");
                par.set_attr("classPath", "jmt.engine.random.ExponentialPar");
                par.set_attr("name", "distrPar");
                jmt_scalar(par, "java.lang.Double", "lambda", "1.000000000000");
                continue;
            }
            jmt_append_distribution(sts, jmt_dist_view(it->second.retrial_proc[r - 1]),
                                    "SolverJMT (retrial delay)");
        }
    }

    /**
     * Port of the polling branch of `saveSwitchoverStrategy`.
     *
     * ONLY the polling branch is ported: `writeJSIM` reaches this handler from
     * the PollingServer section alone, and the reference's second branch --
     * a (K x K) per-class-pair switchover for an ordinary Server -- is
     * unreachable there. A switchover declared on a non-polling queue is
     * refused in `warn_switchover_on_non_polling`, which is where the reference
     * warns and drops it.
     */
    void save_switchover_strategy(xml::Element& section, std::size_t ind) {
        const std::size_t ist = sn_.nodes[ind - 1].station;
        const typename qn::NetworkStruct<T>::PollingParam pp = sn_.effective_polling(ist);
        xml::Element& p = jmt_param(section, "jmt.engine.NetStrategies.ServiceStrategy",
                                    "SwitchoverStrategy", true);
        for (std::size_t r = 1; r <= sn_.nclasses; ++r) {
            jmt_ref_class(p, cname(r));
            xml::Element& sts = p.add_child("subParameter");
            sts.set_attr("classPath",
                         "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            sts.set_attr("name", "ServiceTimeStrategy");
            const lang::Distrib<T>& so =
                r <= pp.switchover.size() ? pp.switchover[r - 1] : empty_dist_;
            const JmtDistView<T> v = jmt_dist_view(so);
            if (so.disabled || v.type == lang::ProcessType::IMMEDIATE) {
                // No switchover declared for this buffer is a zero switchover,
                // not an omission: JMT reads the array positionally. An
                // Immediate switchover is the same zero, declared: the
                // reference maps both to ZeroServiceTimeStrategy, which carries
                // no distribution, and there is no JMT distribution to fall
                // back on -- an unhandled Immediate aborted the whole solve.
                sts.set_attr("classPath",
                             "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                sts.set_attr("name", "ZeroServiceTimeStrategy");
                continue;
            }
            jmt_append_distribution(sts, v, "SolverJMT (switchover)");
        }
    }

    std::map<std::size_t, double> empty_weights_;
    lang::Distrib<T> empty_dist_;
};

/** Port of `@@JMTIO/writeJSIM.m`: serialize `sn` as a JMT `.jsimg` document. */
template <class T>
std::string jmt_write_jsim(const qn::NetworkStruct<T>& sn, const JmtWriteOptions& opt) {
    JmtWriter<T> w(sn, opt);
    return w.write_jsim();
}

/**
 * `JmtWriter::buffer_capacity_refusal` without a document: the writer's own
 * binding-buffer verdict, for the gate in `jmt::jmt_method_refusal`.
 */
template <class T>
std::string jmt_buffer_capacity_refusal(const qn::NetworkStruct<T>& sn, bool jmva_engine) {
    JmtWriter<T> w(sn, JmtWriteOptions());
    return w.buffer_capacity_refusal(jmva_engine);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_JMT_WRITER_H
