/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_JMVA_WRITER_H
#define LINE_IO_JMVA_WRITER_H

/**
 * Port of `@@JMTIO/writeJMVA.m`: the CHAIN-level product-form model in the JMVA
 * interchange format.
 *
 * The format predates this port by way of JMT's MVA panel, and `qnsolver` of the
 * LQNS distribution reads the same grammar, which is why the writer lives in
 * `io/` rather than under either solver: SolverQNS is its only caller here, but
 * the document it writes is the JMT one, field for field.
 *
 * WHAT IT WRITES IS ALREADY AGGREGATED. Classes become CHAINS, service times and
 * visits come from `sn_get_demands_chain`, and a station that is neither a Queue
 * nor a Delay is dropped -- a Source contributes its arrival rate to the open
 * chain's `rate` attribute and nothing else, and the Sink is not a station. The
 * caller therefore gets chain-level results back and has to de-aggregate them;
 * nothing in this file is per-class.
 *
 * A MULTISERVER QUEUE IS WRITTEN AS A LOAD-DEPENDENT STATION, not as a station
 * with `servers` set: the format's own multiserver support is thinner than the
 * rate vector S/min(n,c), and the reference has always spelt the rates out. The
 * `servers` attribute is still emitted as "1" because the schema requires it.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace detail {

/**
 * Shortest decimal text that reads back as the same double.
 *
 * MATLAB writes these through `num2str`, which keeps about five significant
 * digits, and the JAR through `String.valueOf`. Neither is a deliberate choice
 * of precision, and truncating a service demand before handing it to an external
 * solver loses accuracy the port has no way to recover, so the round-trip form
 * is used here.
 */
inline std::string num_text(double x) {
    if (std::isnan(x)) return "0";
    if (x == std::floor(x) && std::fabs(x) < 1e15) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.0f", x);
        return std::string(buf);
    }
    char buf[64];
    for (int prec = 15; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof(buf), "%.*g", prec, x);
        if (std::strtod(buf, nullptr) == x) break;
    }
    return std::string(buf);
}

/** XML character-data and attribute-value escaping. */
inline std::string xml_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        switch (s[i]) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out += s[i];
        }
    }
    return out;
}

/** `sprintf('Chain%02d', c)` of the reference, with c 1-based. */
inline std::string chain_name(std::size_t c1) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "Chain%02zu", c1);
    return std::string(buf);
}

/**
 * The `algType` name the reference maps `options.method` to, and whether that
 * algorithm admits a multiserver station.
 *
 * The `jmva.*` names are JMT's own solvers and reach this writer through
 * SolverJMT; SolverQNS passes a multiserver approximation name, which falls
 * through to plain MVA because the algorithm there is chosen by `qnsolver -m`
 * and not by the document.
 */
inline std::string alg_type_name(const std::string& method, bool* multiserver_ok) {
    *multiserver_ok = true;
    if (method == "jmva.recal") { *multiserver_ok = false; return "RECAL"; }
    if (method == "jmva.comom") { *multiserver_ok = false; return "CoMoM"; }
    if (method == "jmva.chow") { *multiserver_ok = false; return "Chow"; }
    if (method == "jmva.bs" || method == "jmva.amva") {
        *multiserver_ok = false;
        return "Bard-Schweitzer";
    }
    if (method == "jmva.aql") { *multiserver_ok = false; return "AQL"; }
    if (method == "jmva.lin") { *multiserver_ok = false; return "Linearizer"; }
    if (method == "jmva.dmlin") {
        *multiserver_ok = false;
        return "De Souza-Muntz Linearizer";
    }
    return "MVA";
}

}  // namespace detail

/**
 * Port of `writeJMVA(sn, outputFileName, options)`.
 *
 * @param L        the refreshed struct
 * @param path     the file to write
 * @param method   `options.method`, which selects the `algType` name
 * @param samples  `options.samples`, the `maxSamples` attribute
 * @return the path written, so the caller can chain it as the reference does
 */
template <class T>
std::string write_jmva(const qn::NetworkStruct<T>& L, const std::string& path,
                       const std::string& method, std::size_t samples) {
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;

    bool multiserver_ok = true;
    const std::string algname = detail::alg_type_name(method, &multiserver_ok);
    if (!multiserver_ok) {
        for (std::size_t i = 0; i < M; ++i) {
            const double c = L.stations[i].nservers;
            if (std::isfinite(c) && c > 1.0)
                throw UnsupportedError("writeJMVA: " + method +
                                       " does not support multi-server stations");
        }
    }

    const mva::ChainDemands<T> d = mva::sn_get_demands_chain(L);
    auto dv = [](const Matrix<T>& A, std::size_t i, std::size_t j) {
        return num_traits<T>::to_double(A(i, j));
    };

    // The reference indexes sn.rates with a logical over NODES to pick the
    // Source row, which only lands on the right row because the Source is
    // created before any node that is not a station. The station's own node type
    // is the index-safe spelling of the same test and agrees wherever the
    // reference works.
    std::vector<bool> is_source_station(M, false);
    std::size_t nsources = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source) {
            is_source_station[i] = true;
            ++nsources;
        }

    std::ofstream f(path.c_str());
    if (!f) throw InputError("writeJMVA: cannot open '" + path + "' for writing");
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    f << "<model xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
      << " xsi:noNamespaceSchemaLocation=\"JMTmodel.xsd\">\n";
    f << "  <parameters>\n";

    // ---- classes: one per CHAIN ------------------------------------------
    f << "    <classes number=\"" << C << "\">\n";
    for (std::size_t c = 0; c < C; ++c) {
        double sum_njobs = 0.0;
        for (std::size_t k : L.inchain[c]) sum_njobs += L.classes[k - 1].population;
        if (std::isfinite(sum_njobs)) {
            f << "      <closedclass population=\"" << detail::num_text(d.Nchain[c])
              << "\" name=\"" << detail::chain_name(c + 1) << "\"/>\n";
        } else {
            double rate = 0.0;
            for (std::size_t i = 0; i < M; ++i) {
                if (!is_source_station[i]) continue;
                for (std::size_t k : L.inchain[c]) {
                    const double r = dv(L.rates, i, k - 1);
                    if (std::isfinite(r)) rate += r;
                }
            }
            f << "      <openclass rate=\"" << detail::num_text(rate) << "\" name=\""
              << detail::chain_name(c + 1) << "\"/>\n";
        }
    }
    f << "    </classes>\n";

    // ---- stations: Queue and Delay only ----------------------------------
    f << "    <stations number=\"" << (M - nsources) << "\">\n";
    // Whether the closed population is finite decides how far the rate vector of
    // a load-dependent station has to run: an open chain never bounds it, so the
    // vector stops at the server count, past which S/min(n,c) is constant.
    bool any_open = false;
    for (std::size_t k = 0; k < K; ++k)
        if (!std::isfinite(L.classes[k].population)) any_open = true;
    double total_njobs = 0.0;
    for (std::size_t k = 0; k < K; ++k)
        if (std::isfinite(L.classes[k].population)) total_njobs += L.classes[k].population;

    for (std::size_t i = 0; i < M; ++i) {
        const qn::NodeType nt = L.stations[i].nodetype;
        if (nt != qn::NodeType::Queue && nt != qn::NodeType::Delay) continue;
        const std::string name =
            detail::xml_escape(L.nodes[L.station_to_node[i] - 1].name);
        // Effective server count. A load-dependent scaling reaches JMVA as the c
        // of an <ldstation>, the same encoding `save_number_of_servers` uses for
        // JSIM: `check_model` admits only alpha(n) = min(n,c), so max(alpha) is
        // that c. Reading `nservers` alone wrote a <listation> at nominal service
        // time and dropped the scaling.
        double servers = L.stations[i].nservers;
        for (const T& s : L.stations[i].lldscaling)
            servers = std::max(servers, num_traits<T>::to_double(s));
        const bool is_ld = (nt == qn::NodeType::Queue) && !(servers == 1.0);
        const char* tag = nt == qn::NodeType::Delay
                              ? "delaystation"
                              : (is_ld ? "ldstation" : "listation");

        f << "      <" << tag << " name=\"" << name << "\"";
        if (nt == qn::NodeType::Queue) f << " servers=\"1\"";
        f << ">\n";

        f << "        <servicetimes>\n";
        for (std::size_t c = 0; c < C; ++c) {
            const double st = dv(d.STchain, i, c);
            if (is_ld) {
                const double limit = any_open ? servers : total_njobs;
                std::string s = detail::num_text(st);
                for (double n = 2.0; n <= limit; n += 1.0)
                    s += ";" + detail::num_text(st / std::min(n, servers));
                f << "          <servicetimes customerclass=\"" << detail::chain_name(c + 1)
                  << "\">" << s << "</servicetimes>\n";
            } else {
                f << "          <servicetime customerclass=\"" << detail::chain_name(c + 1)
                  << "\">" << detail::num_text(st) << "</servicetime>\n";
            }
        }
        f << "        </servicetimes>\n";

        f << "        <visits>\n";
        for (std::size_t c = 0; c < C; ++c) {
            const double st = dv(d.STchain, i, c);
            const double v = st > 0.0 ? dv(d.Lchain, i, c) / st : 0.0;
            f << "          <visit customerclass=\"" << detail::chain_name(c + 1) << "\">"
              << detail::num_text(v) << "</visit>\n";
        }
        f << "        </visits>\n";
        f << "      </" << tag << ">\n";
    }
    f << "    </stations>\n";

    // ---- reference stations ----------------------------------------------
    // An open chain's reference station is the Source, which is not in the
    // document at all, so the reference substitutes the first station that is
    // neither Source nor Sink. Naming an absent station makes qnsolver reject
    // the whole model, so this is not cosmetic.
    f << "    <ReferenceStation number=\"" << C << "\">\n";
    for (std::size_t c = 0; c < C; ++c) {
        std::size_t ref = L.classes[L.inchain[c][0] - 1].refstat;  // 1-based station
        if (L.stations[ref - 1].nodetype == qn::NodeType::Source) {
            for (std::size_t i = 0; i < M; ++i) {
                const qn::NodeType nt = L.stations[i].nodetype;
                if (nt != qn::NodeType::Source && nt != qn::NodeType::Sink) {
                    ref = i + 1;
                    break;
                }
            }
        }
        f << "      <Class name=\"" << detail::chain_name(c + 1) << "\" refStation=\""
          << detail::xml_escape(L.nodes[L.station_to_node[ref - 1] - 1].name) << "\"/>\n";
    }
    f << "    </ReferenceStation>\n";
    f << "  </parameters>\n";

    f << "  <algParams>\n";
    f << "    <algType name=\"" << detail::xml_escape(algname)
      << "\" tolerance=\"1.0E-7\" maxSamples=\"" << samples << "\"/>\n";
    f << "    <compareAlgs value=\"false\"/>\n";
    f << "  </algParams>\n";
    f << "</model>\n";
    f.close();
    if (!f) throw InputError("writeJMVA: failed to write '" + path + "'");
    return path;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_JMVA_WRITER_H
