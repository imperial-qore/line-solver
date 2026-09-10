/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_PRINT_H
#define LINE_API_SN_SN_PRINT_H

/**
 * Port of matlab/src/api/sn/sn_print.m, jline.api.sn.SnPrint and the python
 * `sn_print` in api/sn/utils.py: the full debug dump of a NetworkStruct, one
 * `field: value` line per field, matrices in the compact `[a b; c d]` form,
 * integer-valued entries printed as integers.
 *
 * The field list is this port's own: the C++ NetworkStruct is object-shaped
 * (per-station and per-class records instead of MATLAB's flat matrices), so
 * the dump renders the MATLAB fields it can derive (refstat, njobs, nservers,
 * phases, chains, ...) in the reference's order and naming, then the
 * structural members that only exist here (station_to_node, stateful_nodes).
 * MATLAB fields with no counterpart in this port (connmatrix, spaceHash) are
 * omitted rather than printed empty; `network_struct.h` documents why each is
 * absent.
 *
 * Returned as a string rather than printed, so a caller can route it to a
 * log, a CLI or a test, following `sn_print_routing_matrix`.
 *
 * ARITHMETIC: field, but the rendering is in double.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

namespace sn_print_detail {

/** One number, integer-rendered when it is one, as the reference prints. */
inline std::string num(double v) {
    if (std::isnan(v)) return "NaN";
    if (std::isinf(v)) return v > 0 ? "Inf" : "-Inf";
    if (v == std::floor(v) && std::fabs(v) < 1e15) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%lld", static_cast<long long>(v));
        return buf;
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%g", v);
    return buf;
}

/** The compact `[a b; c d]` form of the reference's printMatrixCompact. */
inline std::string mat(const std::vector<std::vector<double>>& m) {
    if (m.empty() || m[0].empty()) return "[]";
    std::string out = "[";
    for (std::size_t i = 0; i < m.size(); ++i) {
        if (i > 0) out += "; ";
        for (std::size_t j = 0; j < m[i].size(); ++j) {
            if (j > 0) out += " ";
            out += num(m[i][j]);
        }
    }
    return out + "]";
}

inline std::string row(const std::vector<double>& v) {
    return mat(std::vector<std::vector<double>>(1, v));
}

template <class T>
std::string tmat(const Matrix<T>& m) {
    std::vector<std::vector<double>> d(m.rows(), std::vector<double>(m.cols()));
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) d[i][j] = num_traits<T>::to_double(m(i, j));
    return mat(d);
}

inline std::string names(const std::vector<std::string>& v) {
    if (v.empty()) return "[]";
    std::string out = "[";
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i > 0) out += ", ";
        out += "\"" + v[i] + "\"";
    }
    return out + "]";
}

}  // namespace sn_print_detail

/** The full `field: value` dump of a NetworkStruct, ending with a newline. */
template <class T>
std::string sn_print(const qn::NetworkStruct<T>& sn) {
    using sn_print_detail::mat;
    using sn_print_detail::names;
    using sn_print_detail::num;
    using sn_print_detail::row;
    using sn_print_detail::tmat;
    typedef std::vector<double> Row;
    typedef std::vector<std::vector<double>> Tab;
    const std::size_t M = sn.nof_stations(), K = sn.nclasses, I = sn.nof_nodes();
    std::string out;

    out += "nstations: " + num(static_cast<double>(M)) + "\n";
    out += "nstateful: " + num(static_cast<double>(sn.nof_stateful())) + "\n";
    out += "nnodes: " + num(static_cast<double>(I)) + "\n";
    out += "nclasses: " + num(static_cast<double>(K)) + "\n";
    out += "nclosedjobs: " + num(sn.nclosedjobs()) + "\n";
    out += "nchains: " + num(static_cast<double>(sn.nchains)) + "\n";

    Row refstat(K), njobsv = sn.njobs(), classprio(K);
    for (std::size_t r = 0; r < K; ++r) {
        refstat[r] = static_cast<double>(sn.classes[r].refstat);
        classprio[r] = static_cast<double>(sn.classes[r].prio);
    }
    out += "refstat: " + row(refstat) + "\n";
    out += "njobs: " + row(njobsv) + "\n";
    Row nservers(M), cap(M);
    Tab classcap(M, Row(K)), phases(M, Row(K)), schedparam(M, Row(K, 0.0));
    Tab droprule(M, Row(K)), procid(M, Row(K));
    for (std::size_t i = 0; i < M; ++i) {
        nservers[i] = sn.stations[i].nservers;
        cap[i] = i < sn.cap.size() ? sn.cap[i] : std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < K; ++r) {
            classcap[i][r] = i < sn.classcap.size() && r < sn.classcap[i].size()
                                 ? sn.classcap[i][r]
                                 : std::numeric_limits<double>::infinity();
            phases[i][r] = static_cast<double>(sn.phases_of(i + 1, r + 1));
            if (r < sn.stations[i].schedparam.size())
                schedparam[i][r] = num_traits<T>::to_double(sn.stations[i].schedparam[r]);
            droprule[i][r] = i < sn.droprule.size() && r < sn.droprule[i].size()
                                 ? static_cast<double>(static_cast<int>(sn.droprule[i][r]))
                                 : static_cast<double>(static_cast<int>(qn::DropStrategy::WAITQ));
            procid[i][r] = static_cast<double>(static_cast<int>(sn.procid(i + 1, r + 1)));
        }
    }
    out += "nservers: " + row(nservers) + "\n";
    out += "rates: " + tmat(sn.rates) + "\n";
    out += "scv: " + tmat(sn.scv) + "\n";
    out += "classprio: " + row(classprio) + "\n";
    out += "phases: " + mat(phases) + "\n";
    out += "schedparam: " + mat(schedparam) + "\n";

    Tab chains(sn.chains.size());
    for (std::size_t c = 0; c < sn.chains.size(); ++c)
        for (std::size_t r = 0; r < sn.chains[c].size(); ++r)
            chains[c].push_back(sn.chains[c][r] ? 1.0 : 0.0);
    out += "chains: " + mat(chains) + "\n";
    out += "rt: " + tmat(sn.rt) + "\n";
    out += "rtnodes: " + tmat(sn.rtnodes) + "\n";

    Tab nvars(sn.nvars.size());
    for (std::size_t i = 0; i < sn.nvars.size(); ++i)
        for (std::size_t v = 0; v < sn.nvars[i].size(); ++v)
            nvars[i].push_back(static_cast<double>(sn.nvars[i][v]));
    out += "nvars: " + mat(nvars) + "\n";
    out += "cap: " + row(cap) + "\n";
    out += "classcap: " + mat(classcap) + "\n";

    Row refclass(sn.refclass.size());
    for (std::size_t c = 0; c < sn.refclass.size(); ++c)
        refclass[c] = static_cast<double>(sn.refclass[c]);
    out += "refclass: " + row(refclass) + "\n";

    Tab lld;
    for (std::size_t i = 0; i < M; ++i)
        if (!sn.stations[i].lldscaling.empty()) {
            Row rw;
            for (std::size_t n = 0; n < sn.stations[i].lldscaling.size(); ++n)
                rw.push_back(num_traits<T>::to_double(sn.stations[i].lldscaling[n]));
            lld.push_back(rw);
        }
    out += "lldscaling: " + (lld.empty() ? std::string("[]") : mat(lld)) + "\n";

    Tab fj(sn.fj.size());
    for (std::size_t e = 0; e < sn.fj.size(); ++e) {
        fj[e].push_back(static_cast<double>(sn.fj[e].first));
        fj[e].push_back(static_cast<double>(sn.fj[e].second));
    }
    out += "fj: " + mat(fj) + "\n";

    out += "nodetype: ";
    if (sn.nodes.empty()) {
        out += "[]\n";
    } else {
        out += "[";
        for (std::size_t i = 0; i < I; ++i) {
            if (i > 0) out += ", ";
            out += lang::node_type_to_text(sn.nodes[i].nodetype);
        }
        out += "]\n";
    }
    std::vector<std::string> classnames(K), nodenames(I), schednames(M);
    for (std::size_t r = 0; r < K; ++r) classnames[r] = sn.classes[r].name;
    for (std::size_t i = 0; i < I; ++i) nodenames[i] = sn.nodes[i].name;
    out += "classnames: " + names(classnames) + "\n";
    out += "nodenames: " + names(nodenames) + "\n";

    out += "sched: {";
    for (std::size_t i = 0; i < M; ++i) {
        if (i > 0) out += ", ";
        out += "\"" + sn.stations[i].name + "\": \"" +
               lang::sched_to_text(sn.stations[i].sched) + "\"";
    }
    out += "}\n";

    out += "procid: " + mat(procid) + "\n";
    out += "proc: {";
    for (std::size_t i = 0; i < M; ++i) {
        if (i > 0) out += ", ";
        out += "\"" + sn.stations[i].name + "\": {";
        for (std::size_t r = 0; r < K; ++r) {
            if (r > 0) out += ", ";
            const lang::Distrib<T>& d = sn.service[i][r];
            out += "\"" + sn.classes[r].name + "\": ";
            if (d.disabled) {
                out += "null";
            } else {
                out += std::string("{\"type\": \"") + lang::process_to_text(d.type) +
                       "\", \"mean\": " + num(num_traits<T>::to_double(d.mean)) +
                       ", \"scv\": " + num(num_traits<T>::to_double(d.scv)) + "}";
            }
        }
        out += "}";
    }
    out += "}\n";

    out += "inchain: {";
    for (std::size_t c = 0; c < sn.inchain.size(); ++c) {
        if (c > 0) out += ", ";
        Row rw;
        for (std::size_t k = 0; k < sn.inchain[c].size(); ++k)
            rw.push_back(static_cast<double>(sn.inchain[c][k]));
        out += num(static_cast<double>(c)) + ": " + row(rw);
    }
    out += "}\n";

    out += "visits: {";
    for (std::size_t c = 0; c < sn.visits.size(); ++c) {
        if (c > 0) out += ", ";
        out += num(static_cast<double>(c)) + ": " + tmat(sn.visits[c]);
    }
    out += "}\n";
    out += "nodevisits: {";
    for (std::size_t c = 0; c < sn.nodevisits.size(); ++c) {
        if (c > 0) out += ", ";
        out += num(static_cast<double>(c)) + ": " + tmat(sn.nodevisits[c]);
    }
    out += "}\n";
    out += "droprule: " + mat(droprule) + "\n";

    Row s2n(sn.station_to_node.size()), stf(sn.stateful_nodes.size());
    for (std::size_t i = 0; i < sn.station_to_node.size(); ++i)
        s2n[i] = static_cast<double>(sn.station_to_node[i]);
    for (std::size_t i = 0; i < sn.stateful_nodes.size(); ++i)
        stf[i] = static_cast<double>(sn.stateful_nodes[i]);
    out += "stationToNode: " + row(s2n) + "\n";
    out += "statefulNodes: " + row(stf) + "\n";

    out += "csmatrix: {";
    {
        bool first = true;
        for (typename std::map<std::size_t, Matrix<T>>::const_iterator it = sn.csmatrix.begin();
             it != sn.csmatrix.end(); ++it) {
            if (!first) out += ", ";
            first = false;
            out += num(static_cast<double>(it->first)) + ": " + tmat(it->second);
        }
    }
    out += "}\n";
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_PRINT_H
