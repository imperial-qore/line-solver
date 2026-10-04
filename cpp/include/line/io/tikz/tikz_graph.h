/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_TIKZ_TIKZ_GRAPH_H
#define LINE_IO_TIKZ_TIKZ_GRAPH_H

/**
 * @file
 * @ingroup line_io
 * The part of a network the TikZ exporter draws, and the text helpers the
 * JAR's exporter uses to write it.
 *
 * `jline.io.tikz` reads a live `Network`: its node objects (for their class and
 * name), each Queue's discipline and server count, and `sn.connmatrix`. This
 * port reads the same facts off a `NetworkStruct` once, into a `TikzGraph`, so
 * the layout, routing and rendering stages below are plain functions of it.
 *
 * THE CONNECTION MATRIX IS THE UNION SUPPORT OF THE ROUTING BLOCKS `P`, as in
 * `jmt_conn_matrix` (jmt_writer.h): this port records links only through `P`,
 * and `Network::link` has already replaced every class-switching link by the
 * `CS_i_to_j` node the JAR inserts, so the topology is the JAR's, node for node.
 *
 * NUMBERS ARE FORMATTED AS JAVA'S `String.format("%.2f")` DOES, not as printf
 * does. Java rounds HALF UP on the shortest decimal expansion of the double,
 * printf rounds the exact binary value half to even, and the two disagree on a
 * coordinate such as 1.675 (printf 1.67, Java 1.68). Such coordinates do occur:
 * a forward edge is routed at the MEAN height of its obstacles.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace io {

/** One node as the TikZ exporter sees it. */
struct TikzNode {
    std::string name;
    lang::NodeType type = lang::NodeType::Queue;
    std::string sched;     ///< Java `SchedStrategy.name()` of a Queue (FCFSPRIO is HOL, as in MATLAB), empty when none
    double servers = 1.0;  ///< a Queue's server count, infinite for an unbounded pool
};

/** Nodes in model order plus `sn.connmatrix` over them. */
struct TikzGraph {
    std::string name;
    std::vector<TikzNode> nodes;
    std::vector<std::vector<bool>> conn;  ///< conn[i][j]: node i is linked to node j, 0-based
};

namespace tikz_detail {

/** Java's `SchedStrategy.name()`: the enumerator spelling, which is the upper-cased text of `sched_to_text`. */
inline std::string java_sched_name(lang::SchedStrategy s) {
    if (s == lang::SchedStrategy::NONE) return std::string();
    std::string t = lang::sched_to_text(s);
    for (std::size_t k = 0; k < t.size(); ++k)
        if (t[k] >= 'a' && t[k] <= 'z') t[k] = static_cast<char>(t[k] - 'a' + 'A');
    return t;
}

/**
 * Java's `String.format("%." + prec + "f", v)`.
 *
 * The digits are the shortest decimal expansion that reads back as `v` (what
 * `Double.toString` prints), rounded half up at `prec` decimals, which is
 * `FormattedFloatingDecimal.applyPrecision`. A negative value keeps its sign
 * even when it rounds to zero, as Java's does ("-0.00").
 */
inline std::string java_fixed(double v, int prec) {
    if (std::isnan(v)) return "NaN";
    const bool neg = std::signbit(v);
    const double a = std::fabs(v);
    if (std::isinf(a)) return neg ? "-Infinity" : "Infinity";
    std::string digits;  // significant digits, value = 0.d1d2... * 10^pt
    int pt = 1;
    if (a == 0.0) {
        digits = "0";
    } else {
        char buf[64];
        for (int p = 1; p <= 17; ++p) {
            std::snprintf(buf, sizeof(buf), "%.*e", p - 1, a);
            if (std::strtod(buf, nullptr) == a) break;
        }
        const std::string s(buf);
        const std::size_t e = s.find('e');
        for (std::size_t k = 0; k < e; ++k)
            if (s[k] != '.') digits.push_back(s[k]);
        pt = std::atoi(s.c_str() + e + 1) + 1;
    }
    if (pt <= 0) {
        digits = std::string(static_cast<std::size_t>(1 - pt), '0') + digits;
        pt = 1;
    }
    const std::size_t keep = static_cast<std::size_t>(pt + prec);
    const bool round_up = digits.size() > keep && digits[keep] >= '5';
    if (digits.size() < keep) digits.append(keep - digits.size(), '0');
    digits.resize(keep);
    if (round_up) {
        std::size_t k = keep;
        while (k > 0) {
            --k;
            if (digits[k] == '9') {
                digits[k] = '0';
            } else {
                ++digits[k];
                break;
            }
            if (k == 0) {
                digits.insert(digits.begin(), '1');
                ++pt;
            }
        }
    }
    std::string ip = digits.substr(0, static_cast<std::size_t>(pt));
    const std::size_t nz = ip.find_first_not_of('0');
    ip = (nz == std::string::npos) ? std::string("0") : ip.substr(nz);
    std::string out = neg ? "-" : "";
    out += ip;
    if (prec > 0) out += "." + digits.substr(static_cast<std::size_t>(pt));
    return out;
}

/** `name.replaceAll("[^a-zA-Z0-9]", "_")`: one underscore per CODE POINT, so a UTF-8 sequence is one character. */
inline std::string sanitize_id(const std::string& name) {
    std::string out;
    out.reserve(name.size());
    for (std::size_t k = 0; k < name.size(); ++k) {
        const unsigned char c = static_cast<unsigned char>(name[k]);
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9')) {
            out.push_back(static_cast<char>(c));
        } else if ((c & 0xC0) != 0x80) {
            out.push_back('_');  // a lead byte or any other ASCII character; continuation bytes add nothing
        }
    }
    return out;
}

inline void replace_all(std::string& s, const std::string& from, const std::string& to) {
    std::size_t at = 0;
    while ((at = s.find(from, at)) != std::string::npos) {
        s.replace(at, from.size(), to);
        at += to.size();
    }
}

/**
 * `TikZNodeRenderer.escapeLatex`, applied in the SAME sequence of whole-string
 * replacements, so a backslash comes out as `\textbackslash\{\}`, exactly as
 * the JAR's later brace replacement leaves it.
 */
inline std::string escape_latex(std::string t) {
    replace_all(t, "\\", "\\textbackslash{}");
    replace_all(t, "_", "\\_");
    replace_all(t, "&", "\\&");
    replace_all(t, "%", "\\%");
    replace_all(t, "$", "\\$");
    replace_all(t, "#", "\\#");
    replace_all(t, "{", "\\{");
    replace_all(t, "}", "\\}");
    replace_all(t, "~", "\\textasciitilde{}");
    replace_all(t, "^", "\\textasciicircum{}");
    return t;
}

}  // namespace tikz_detail

/**
 * The drawable graph of a network struct.
 *
 * A node is drawn by the CLASS it was declared with, as MATLAB (through the JAR's
 * `instanceof` chain) draws it: a Queue scheduled INF has nodetype Delay in the
 * struct (`Network::add_queue`), but `queue_object` marks it, so it is drawn as
 * a buffer labelled INF with an infinite server, not as a Delay box.
 */
template <class T>
TikzGraph tikz_graph(const qn::NetworkStruct<T>& sn) {
    TikzGraph g;
    g.name = sn.name;
    const std::size_t I = sn.nodes.size();
    g.nodes.resize(I);
    for (std::size_t i = 0; i < I; ++i) {
        const qn::NodeDef& nd = sn.nodes[i];
        TikzNode& tn = g.nodes[i];
        tn.name = nd.name;
        tn.type = nd.queue_object ? lang::NodeType::Queue : nd.nodetype;
        if (nd.station > 0 && nd.station <= sn.stations.size()) {
            const qn::Station<T>& st = sn.stations[nd.station - 1];
            tn.sched = tikz_detail::java_sched_name(st.sched);
            tn.servers = st.nservers;
        }
    }
    g.conn.assign(I, std::vector<bool>(I, false));
    const T zero = num_traits<T>::from_int(0);
    for (const auto& kv : sn.P) {
        const Matrix<T>& B = kv.second;
        for (std::size_t a = 0; a < B.rows() && a < I; ++a)
            for (std::size_t b = 0; b < B.cols() && b < I; ++b)
                if (!(B(a, b) == zero)) g.conn[a][b] = true;
    }
    // The Sink -> Source arc is the refresh's closure of the open chains (`apply_sink_closure`), which the
    // reference builds on the routing matrix and never records as a link; it is not part of the drawing.
    if (sn.sourceIdx > 0 && sn.sinkNode > 0 && sn.sourceIdx <= sn.station_to_node.size()) {
        const std::size_t src = sn.station_to_node[sn.sourceIdx - 1];
        if (sn.sinkNode <= I && src >= 1 && src <= I) g.conn[sn.sinkNode - 1][src - 1] = false;
    }
    return g;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_TIKZ_TIKZ_GRAPH_H
