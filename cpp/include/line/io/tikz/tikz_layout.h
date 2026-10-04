/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_TIKZ_TIKZ_LAYOUT_H
#define LINE_IO_TIKZ_TIKZ_LAYOUT_H

/**
 * @file
 * @ingroup line_io
 * Port of `jline.io.tikz.TikZLayoutEngine`: a layered (Sugiyama-style) layout.
 *
 * Three steps, each the JAR's: layers by a breadth-first topological sweep from
 * the Sources and the nodes nobody feeds (a cycle is broken by taking the first
 * unassigned node in model order), Sinks moved to the last layer; four
 * forward/backward barycenter passes to reduce crossings, each a STABLE sort as
 * `Collections.sort` is; then x = layer * layerSpacing and each layer centred
 * on y = 0 with nodeSpacing between its members.
 */

#include <algorithm>
#include <limits>
#include <vector>

#include "line/io/tikz/tikz_graph.h"
#include "line/io/tikz/tikz_options.h"

namespace line {
namespace io {

/** Node positions in cm, indexed like `TikzGraph::nodes`, and the layers that produced them. */
struct TikzLayout {
    std::vector<double> x, y;
    std::vector<std::vector<std::size_t>> layers;  ///< 0-based node indices, in drawing order
};

namespace tikz_detail {

/** `reorderLayerByBarycenter`: sort a layer by the mean position of its neighbours in the reference layer. */
inline void reorder_by_barycenter(std::vector<std::size_t>& layer,
                                  const std::vector<std::vector<std::size_t>>& connections,
                                  const std::vector<std::size_t>& reference, std::size_t nnodes) {
    std::vector<long> refpos(nnodes, -1);
    for (std::size_t i = 0; i < reference.size(); ++i) refpos[reference[i]] = static_cast<long>(i);
    std::vector<double> bary(nnodes, std::numeric_limits<double>::max());
    for (std::size_t k = 0; k < layer.size(); ++k) {
        const std::size_t v = layer[k];
        double sum = 0;
        int count = 0;
        for (std::size_t c : connections[v]) {
            if (refpos[c] >= 0) {
                sum += static_cast<double>(refpos[c]);
                ++count;
            }
        }
        if (count > 0) bary[v] = sum / count;
    }
    std::stable_sort(layer.begin(), layer.end(),
                     [&](std::size_t a, std::size_t b) { return bary[a] < bary[b]; });
}

}  // namespace tikz_detail

/** `TikZLayoutEngine.computeLayout`. */
inline TikzLayout tikz_layout(const TikzGraph& g, const TikzOptions& opt) {
    TikzLayout L;
    const std::size_t n = g.nodes.size();
    L.x.assign(n, 0.0);
    L.y.assign(n, 0.0);
    if (n == 0) return L;

    std::vector<std::vector<std::size_t>> succ(n), pred(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (g.conn[i][j]) {
                succ[i].push_back(j);
                pred[j].push_back(i);
            }

    // Step 1: layers
    std::vector<std::vector<std::size_t>>& layers = L.layers;
    std::vector<bool> assigned(n, false);
    std::size_t nassigned = 0;
    std::vector<std::size_t> layer0;
    for (std::size_t v = 0; v < n; ++v)
        if (g.nodes[v].type == lang::NodeType::Source || pred[v].empty()) {
            layer0.push_back(v);
            assigned[v] = true;
            ++nassigned;
        }
    if (layer0.empty()) {
        layer0.push_back(0);
        assigned[0] = true;
        ++nassigned;
    }
    layers.push_back(layer0);
    std::size_t cur = 0;
    while (nassigned < n) {
        std::vector<std::size_t> next;
        for (std::size_t v = 0; v < n; ++v) {
            if (assigned[v]) continue;
            bool all = true;
            for (std::size_t p : pred[v])
                if (!assigned[p]) {
                    all = false;
                    break;
                }
            if (all) next.push_back(v);
        }
        if (next.empty())  // a cycle: take the first unassigned node
            for (std::size_t v = 0; v < n; ++v)
                if (!assigned[v]) {
                    next.push_back(v);
                    break;
                }
        if (!next.empty()) {
            for (std::size_t v : next) {
                assigned[v] = true;
                ++nassigned;
            }
            layers.push_back(next);
        }
        ++cur;
        if (cur > n) break;  // the JAR's own guard
    }
    if (layers.size() > 1) {
        std::vector<std::size_t>& last = layers.back();
        for (std::size_t i = 0; i + 1 < layers.size(); ++i) {
            std::vector<std::size_t>& layer = layers[i];
            for (std::size_t k = 0; k < layer.size();) {
                const std::size_t v = layer[k];
                if (g.nodes[v].type == lang::NodeType::Sink) {
                    layer.erase(layer.begin() + static_cast<long>(k));
                    if (std::find(last.begin(), last.end(), v) == last.end()) last.push_back(v);
                } else {
                    ++k;
                }
            }
        }
        std::vector<std::vector<std::size_t>> kept;
        for (std::size_t i = 0; i < layers.size(); ++i)
            if (!layers[i].empty()) kept.push_back(layers[i]);
        layers.swap(kept);
    }

    // Step 2: crossing reduction
    for (int pass = 0; pass < 4; ++pass) {
        for (std::size_t i = 1; i < layers.size(); ++i)
            tikz_detail::reorder_by_barycenter(layers[i], pred, layers[i - 1], n);
        for (std::size_t i = layers.size() >= 2 ? layers.size() - 1 : 0; i-- > 0;)
            tikz_detail::reorder_by_barycenter(layers[i], succ, layers[i + 1], n);
    }

    // Step 3: coordinates
    for (std::size_t li = 0; li < layers.size(); ++li) {
        const std::vector<std::size_t>& layer = layers[li];
        const double x = static_cast<double>(li) * opt.layer_spacing;
        const double total = static_cast<double>(layer.size() - 1) * opt.node_spacing;
        const double start = total / 2.0;
        for (std::size_t k = 0; k < layer.size(); ++k) {
            L.x[layer[k]] = x;
            L.y[layer[k]] = start - static_cast<double>(k) * opt.node_spacing;
        }
    }
    return L;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_TIKZ_TIKZ_LAYOUT_H
