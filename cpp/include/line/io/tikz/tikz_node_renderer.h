/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_TIKZ_TIKZ_NODE_RENDERER_H
#define LINE_IO_TIKZ_TIKZ_NODE_RENDERER_H

/**
 * @file
 * @ingroup line_io
 * Port of `jline.io.tikz.TikZNodeRenderer`: the preamble with the node styles,
 * and one TikZ node per network node.
 *
 * The shape follows the node's KIND, in the JAR's `instanceof` order: a Delay
 * is a green box marked infinity, a Queue a buffer with its discipline beneath
 * and a server circle carrying the server count, a Source/Sink a small
 * white/black circle with an arrival/departure stub, Fork/Join diamonds, a
 * Router a hexagon, a ClassSwitch a trapezium, a Cache a grey box with three
 * level lines, a Logger a brown box, a Place a circle, a Transition a black bar.
 * Any other kind (a Region) falls back to a plain circle.
 */

#include <cmath>
#include <sstream>
#include <string>

#include "line/io/tikz/tikz_graph.h"
#include "line/io/tikz/tikz_options.h"

namespace line {
namespace io {

/** `TikZNodeRenderer.getPreamble`, byte for byte. */
inline std::string tikz_preamble(const TikzOptions& opt) {
    std::ostringstream sb;
    sb << "\\documentclass[tikz,border=" << opt.border_padding << "pt]{standalone}\n";
    sb << "\\usepackage{tikz}\n";
    sb << "\\usetikzlibrary{arrows.meta,positioning,shapes.geometric,shapes.misc,calc,"
          "decorations.pathreplacing}\n";
    sb << "\n";
    sb << "\\tikzset{\n";
    sb << "    % Queue: Rectangle buffer\n"
          "    queue/.style={\n"
          "        rectangle,\n"
          "        draw=black,\n"
          "        minimum width=1.1cm,\n"
          "        minimum height=0.72cm,\n"
          "        fill=white\n"
          "    },\n";
    sb << "    % Server circle\n"
          "    server/.style={\n"
          "        circle,\n"
          "        draw=black,\n"
          "        minimum size=0.6cm,\n"
          "        fill=white\n"
          "    },\n";
    sb << "    % Delay: Vertical rectangle (infinite server)\n"
          "    delay/.style={\n"
          "        rectangle,\n"
          "        draw=black,\n"
          "        minimum width=0.8cm,\n"
          "        minimum height=1.5cm,\n"
          "        fill=green!10\n"
          "    },\n";
    sb << "    % Source: Small circle (white)\n"
          "    source/.style={\n"
          "        circle,\n"
          "        draw=black,\n"
          "        minimum size=0.3cm,\n"
          "        fill=white\n"
          "    },\n";
    sb << "    % Sink: Small circle (black)\n"
          "    sink/.style={\n"
          "        circle,\n"
          "        draw=black,\n"
          "        minimum size=0.3cm,\n"
          "        fill=black\n"
          "    },\n";
    sb << "    % Fork: Diamond\n"
          "    fork/.style={\n"
          "        diamond,\n"
          "        draw=black,\n"
          "        minimum size=1cm,\n"
          "        fill=orange!20,\n"
          "        aspect=1.5\n"
          "    },\n";
    sb << "    % Join: Diamond\n"
          "    joinnode/.style={\n"
          "        diamond,\n"
          "        draw=black,\n"
          "        minimum size=1cm,\n"
          "        fill=purple!20,\n"
          "        aspect=1.5\n"
          "    },\n";
    sb << "    % Router: Hexagon\n"
          "    router/.style={\n"
          "        regular polygon,\n"
          "        regular polygon sides=6,\n"
          "        draw=black,\n"
          "        minimum size=1cm,\n"
          "        fill=cyan!10\n"
          "    },\n";
    sb << "    % ClassSwitch: Trapezium\n"
          "    classswitch/.style={\n"
          "        trapezium,\n"
          "        draw=black,\n"
          "        trapezium left angle=70,\n"
          "        trapezium right angle=110,\n"
          "        minimum width=1.5cm,\n"
          "        minimum height=0.8cm,\n"
          "        fill=pink!20\n"
          "    },\n";
    sb << "    % Cache: Stacked rectangle\n"
          "    cache/.style={\n"
          "        rectangle,\n"
          "        draw=black,\n"
          "        minimum width=1.5cm,\n"
          "        minimum height=1.2cm,\n"
          "        fill=gray!20\n"
          "    },\n";
    sb << "    % Logger: Rectangle with lines\n"
          "    logger/.style={\n"
          "        rectangle,\n"
          "        draw=black,\n"
          "        minimum width=1.2cm,\n"
          "        minimum height=0.8cm,\n"
          "        fill=brown!10\n"
          "    },\n";
    sb << "    % Place (Petri net): Circle\n"
          "    place/.style={\n"
          "        circle,\n"
          "        draw=black,\n"
          "        minimum size=0.8cm,\n"
          "        fill=white\n"
          "    },\n";
    sb << "    % Transition (Petri net): Rectangle\n"
          "    transition/.style={\n"
          "        rectangle,\n"
          "        draw=black,\n"
          "        minimum width=0.2cm,\n"
          "        minimum height=1cm,\n"
          "        fill=black\n"
          "    },\n";
    sb << "    % Connection arrow\n"
          "    conn/.style={\n"
          "        ->,\n"
          "        >=Stealth,\n"
          "        thick\n"
          "    },\n";
    sb << "    % Probability label\n"
          "    problabel/.style={\n"
          "        font=\\footnotesize,\n"
          "        midway,\n"
          "        above,\n"
          "        sloped\n"
          "    },\n";
    sb << "    % Node name label\n"
          "    nodelabel/.style={\n"
          "        font=\\small\n"
          "    }\n";
    sb << "}\n";
    return sb.str();
}

namespace tikz_detail {

/** `\node[style] (id) at (x,y) {body};` followed by the name label the options ask for. */
inline std::string styled_node(const std::string& style, const std::string& id, double x, double y,
                               const std::string& body, const TikzNode& nd, const TikzOptions& opt) {
    std::string s = "\\node[" + style + "] (" + id + ") at (" + java_fixed(x, 2) + "," +
                    java_fixed(y, 2) + ") {" + body + "};\n";
    if (opt.show_node_names)
        s += "\\node[nodelabel,above=2pt of " + id + "] {" + escape_latex(nd.name) + "};\n";
    return s;
}

/** `String.valueOf(queue.getNumberOfServers())`, the JAR's `Integer.MAX_VALUE` pool being infinity. */
inline std::string server_label(double c) {
    if (std::isinf(c) || c >= 2147483647.0) return "$\\infty$";
    std::ostringstream o;
    o << static_cast<long long>(c);
    return o.str();
}

}  // namespace tikz_detail

/** `TikZNodeRenderer.renderNode`. */
inline std::string tikz_render_node(const TikzNode& nd, double x, double y, const TikzOptions& opt) {
    using lang::NodeType;
    using tikz_detail::styled_node;
    const std::string id = tikz_detail::sanitize_id(nd.name);
    std::string s;
    switch (nd.type) {
        case NodeType::Delay:
            return styled_node("delay", id, x, y, "$\\infty$", nd, opt);
        case NodeType::Queue:
            s = styled_node("queue", id, x, y, "", nd, opt);
            if (opt.show_scheduling && !nd.sched.empty())
                s += "\\node[font=\\tiny,below=2pt of " + id + "] {" + nd.sched + "};\n";
            if (opt.show_server_count)
                s += "\\node[server,anchor=west] (" + id + "_server) at (" + id + ".east) {" +
                     tikz_detail::server_label(nd.servers) + "};\n";
            return s;
        case NodeType::Source:
            s = styled_node("source", id, x, y, "", nd, opt);
            s += "\\draw[conn] ([xshift=-0.6cm]" + id + ".west) -- (" + id + ".west);\n";
            return s;
        case NodeType::Sink:
            s = styled_node("sink", id, x, y, "", nd, opt);
            s += "\\draw[conn] (" + id + ".east) -- ([xshift=0.6cm]" + id + ".east);\n";
            return s;
        case NodeType::Fork:
            return styled_node("fork", id, x, y, "", nd, opt);
        case NodeType::Join:
            return styled_node("joinnode", id, x, y, "", nd, opt);
        case NodeType::Router:
            return styled_node("router", id, x, y, "", nd, opt);
        case NodeType::ClassSwitch:
            return styled_node("classswitch", id, x, y, "", nd, opt);
        case NodeType::Cache:
            s = styled_node("cache", id, x, y, "", nd, opt);
            s += "\\draw ([yshift=-0.3cm]" + id + ".north west) -- ([yshift=-0.3cm]" + id +
                 ".north east);\n";
            s += "\\draw (" + id + ".west) -- (" + id + ".east);\n";
            s += "\\draw ([yshift=0.3cm]" + id + ".south west) -- ([yshift=0.3cm]" + id +
                 ".south east);\n";
            return s;
        case NodeType::Logger:
            return styled_node("logger", id, x, y, "", nd, opt);
        case NodeType::Place:
            return styled_node("place", id, x, y, "", nd, opt);
        case NodeType::Transition:
            return styled_node("transition", id, x, y, "", nd, opt);
        default:
            return styled_node("draw,circle,minimum size=0.8cm", id, x, y, "", nd, opt);
    }
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_TIKZ_TIKZ_NODE_RENDERER_H
