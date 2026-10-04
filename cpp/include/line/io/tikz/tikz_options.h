/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_TIKZ_TIKZ_OPTIONS_H
#define LINE_IO_TIKZ_TIKZ_OPTIONS_H

/**
 * @file
 * @ingroup line_io
 * Port of `jline.io.tikz.TikZOptions`: the knobs of the network TikZ exporter.
 *
 * The defaults are the JAR's constructor defaults, so a default-constructed
 * `TikzOptions` produces the same document as `Network.toTikZ()` with no
 * argument, which is what MATLAB's `MNetwork.toTikZ` calls.
 */

namespace line {
namespace io {

/** Layout and rendering options of the network TikZ exporter (`TikZOptions`). */
struct TikzOptions {
    double node_spacing = 3.0;            ///< vertical distance between nodes of one layer, cm
    double layer_spacing = 4.0;           ///< horizontal distance between layers, cm
    bool show_routing_prob = true;        ///< label edges with their probability (see tikz.h: never reached by the exporter)
    bool show_server_count = true;        ///< draw the server circle of a Queue with its server count
    bool show_node_names = true;          ///< label each node with its name
    bool show_scheduling = true;          ///< print a Queue's discipline beneath it
    int border_padding = 10;              ///< `standalone` border, pt
    double min_prob_to_show = 0.001;      ///< probabilities within this of 0 or 1 get no label
    bool hide_auto_generated_nodes = false;  ///< hide the CS_x_to_y class switches `link` inserts
};

}  // namespace io
}  // namespace line

#endif  // LINE_IO_TIKZ_TIKZ_OPTIONS_H
