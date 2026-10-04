/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_TIKZ_TIKZ_EDGE_ROUTER_H
#define LINE_IO_TIKZ_TIKZ_EDGE_ROUTER_H

/**
 * @file
 * @ingroup line_io
 * Port of `jline.io.tikz.TikZEdgeRouter`: orthogonal waypoints that keep an
 * edge off the nodes it would otherwise cross.
 *
 * A forward edge is straight unless a node's box (plus margin) lies on it, and
 * is then routed above or below the obstacles' mean height. A BACKWARD edge
 * (feedback) is routed around the whole drawing, above it when its source is
 * not below its target, and each one takes the next channel out, so the
 * waypoints of a backward edge depend on how many were routed before it: the
 * router is STATEFUL and the order edges are routed in is part of the output.
 * An edge between two nodes of one layer is offset to the right.
 *
 * The JAR's bounds start at `Double.MIN_VALUE` for the maxima (the smallest
 * POSITIVE double, not the most negative); that is kept, since every drawing
 * has a node at y >= 0 and the two starts then give the same bounds anyway.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/io/tikz/tikz_graph.h"
#include "line/io/tikz/tikz_options.h"

namespace line {
namespace io {

/** A waypoint in cm. */
struct TikzPoint {
    double x = 0, y = 0;
};

class TikzEdgeRouter {
public:
    static constexpr double NODE_HALF_WIDTH = 1.2;   ///< cm
    static constexpr double NODE_HALF_HEIGHT = 0.9;  ///< cm
    static constexpr double MARGIN = 0.25;
    static constexpr double ROUTE_SPACING = 0.5;
    static constexpr double ENTRY_OFFSET = 0.2;

    /**
     * @param x,y     positions of every node of the graph, indexed by node
     * @param visible the nodes that are drawn; only they bound the drawing and count as obstacles
     */
    TikzEdgeRouter(const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<std::size_t>& visible)
        : x_(x), y_(y), visible_(visible) {
        min_x_ = std::numeric_limits<double>::max();
        max_x_ = std::numeric_limits<double>::denorm_min();
        min_y_ = std::numeric_limits<double>::max();
        max_y_ = std::numeric_limits<double>::denorm_min();
        for (std::size_t v : visible_) {
            min_x_ = std::min(min_x_, x_[v] - NODE_HALF_WIDTH);
            max_x_ = std::max(max_x_, x_[v] + NODE_HALF_WIDTH);
            min_y_ = std::min(min_y_, y_[v] - NODE_HALF_HEIGHT);
            max_y_ = std::max(max_y_, y_[v] + NODE_HALF_HEIGHT);
        }
    }

    /** `computeWaypoints`: empty when the straight segment is clear, and for a self-loop. */
    std::vector<TikzPoint> compute_waypoints(std::size_t from, std::size_t to) {
        std::vector<TikzPoint> wp;
        if (from == to) return wp;
        const TikzPoint f{x_[from], y_[from]}, t{x_[to], y_[to]};
        if (t.x < f.x - NODE_HALF_WIDTH) return route_backward(f, t);
        if (std::fabs(t.x - f.x) < NODE_HALF_WIDTH) {
            const double off = NODE_HALF_WIDTH + MARGIN + 0.5;
            wp.push_back(TikzPoint{f.x + off, f.y});
            wp.push_back(TikzPoint{t.x + off, t.y});
            return wp;
        }
        std::vector<std::size_t> obstacles;
        for (std::size_t v : visible_) {
            if (v == from || v == to) continue;
            if (line_intersects_node(f, t, v)) obstacles.push_back(v);
        }
        std::stable_sort(obstacles.begin(), obstacles.end(),
                         [&](std::size_t a, std::size_t b) { return x_[a] < x_[b]; });
        if (obstacles.empty()) return wp;
        return route_forward(f, t, obstacles);
    }

private:
    std::vector<TikzPoint> route_backward(const TikzPoint& f, const TikzPoint& t) {
        std::vector<TikzPoint> wp;
        ++backward_;
        const bool above = f.y >= t.y;
        const int ch = backward_;
        const double route_y = above ? max_y_ + MARGIN + (ROUTE_SPACING * ch)
                                     : min_y_ - MARGIN - (ROUTE_SPACING * ch);
        const double exit_x = f.x + NODE_HALF_WIDTH + MARGIN + (backward_ - 1) * 0.3;
        const double entry_x = t.x - NODE_HALF_WIDTH - MARGIN - (backward_ - 1) * 0.3;
        double dy = (backward_ - 1) * ENTRY_OFFSET;
        if (!above) dy = -dy;
        const double target_y = t.y + dy;
        wp.push_back(TikzPoint{exit_x, f.y});
        wp.push_back(TikzPoint{exit_x, route_y});
        wp.push_back(TikzPoint{entry_x, route_y});
        wp.push_back(TikzPoint{entry_x, target_y});
        return wp;
    }

    std::vector<TikzPoint> route_forward(const TikzPoint& f, const TikzPoint& t,
                                         const std::vector<std::size_t>& obstacles) const {
        std::vector<TikzPoint> wp;
        double avg = 0;
        for (std::size_t v : obstacles) avg += y_[v];
        avg /= static_cast<double>(obstacles.size());
        const double mid = (f.y + t.y) / 2.0;
        const bool above = mid >= avg;
        const double offset = (NODE_HALF_HEIGHT + MARGIN + 0.15) * (above ? 1 : -1);
        const double route_y = avg + offset;
        double minx = std::numeric_limits<double>::max();
        double maxx = std::numeric_limits<double>::denorm_min();
        for (std::size_t v : obstacles) {
            minx = std::min(minx, x_[v] - NODE_HALF_WIDTH - MARGIN);
            maxx = std::max(maxx, x_[v] + NODE_HALF_WIDTH + MARGIN);
        }
        const double entry_x = std::max(f.x + 0.3, minx - 0.3);
        wp.push_back(TikzPoint{entry_x, route_y});
        const double exit_x = std::min(t.x - 0.3, maxx + 0.3);
        if (exit_x > entry_x + 0.1) wp.push_back(TikzPoint{exit_x, route_y});
        return wp;
    }

    bool line_intersects_node(const TikzPoint& s, const TikzPoint& e, std::size_t v) const {
        const double hw = NODE_HALF_WIDTH + MARGIN, hh = NODE_HALF_HEIGHT + MARGIN;
        const double left = x_[v] - hw, right = x_[v] + hw;
        const double bottom = y_[v] - hh, top = y_[v] + hh;
        const double x1 = s.x, y1 = s.y, x2 = e.x, y2 = e.y;
        if ((x1 < left && x2 < left) || (x1 > right && x2 > right)) return false;
        if ((y1 < bottom && y2 < bottom) || (y1 > top && y2 > top)) return false;
        if (in_box(x1, y1, left, right, bottom, top) || in_box(x2, y2, left, right, bottom, top))
            return true;
        return seg(x1, y1, x2, y2, left, bottom, left, top) ||
               seg(x1, y1, x2, y2, right, bottom, right, top) ||
               seg(x1, y1, x2, y2, left, bottom, right, bottom) ||
               seg(x1, y1, x2, y2, left, top, right, top);
    }

    static bool in_box(double x, double y, double l, double r, double b, double t) {
        return x >= l && x <= r && y >= b && y <= t;
    }

    static bool seg(double x1, double y1, double x2, double y2, double x3, double y3, double x4,
                    double y4) {
        const double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (std::fabs(den) < 1e-10) return false;
        const double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
        const double u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / den;
        return t >= 0 && t <= 1 && u >= 0 && u <= 1;
    }

    const std::vector<double>& x_;
    const std::vector<double>& y_;
    std::vector<std::size_t> visible_;
    double min_x_, max_x_, min_y_, max_y_;
    int backward_ = 0;
};

/**
 * `TikZEdgeRouter.renderRoutedEdge`. `prob` is NaN for an unlabelled edge, which
 * is every edge the exporter draws: the JAR passes NaN unconditionally.
 */
inline std::string tikz_render_routed_edge(const std::string& from_id, const std::string& to_id,
                                           const std::string& from_anchor,
                                           const std::string& to_anchor,
                                           const std::vector<TikzPoint>& wp, double prob,
                                           const TikzOptions& opt) {
    using tikz_detail::java_fixed;
    const bool label = opt.show_routing_prob && !std::isnan(prob) &&
                       prob >= opt.min_prob_to_show && prob < 1.0 - opt.min_prob_to_show;
    std::string s = "\\draw[conn] (" + from_id + from_anchor + ")";
    for (const TikzPoint& p : wp) s += " -- (" + java_fixed(p.x, 2) + "," + java_fixed(p.y, 2) + ")";
    if (label) s += " -- node[problabel] {" + java_fixed(prob, 2) + "}";
    else s += " --";
    s += " (" + to_id + to_anchor + ");\n";
    return s;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_TIKZ_TIKZ_EDGE_ROUTER_H
