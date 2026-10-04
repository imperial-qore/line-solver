"""
Orthogonal edge routing around nodes for queueing network diagrams.

Port of jar/src/main/java/jline/io/tikz/TikZEdgeRouter.java. Forward edges are
routed around intermediate nodes, backward (feedback) edges around the whole
network in their own channel, vertical edges with a horizontal offset; self
loops are drawn by the exporter.
"""

import math

from ._jcompat import JAVA_DOUBLE_MAX_VALUE, JAVA_DOUBLE_MIN_VALUE, f2, jmax, jmin

NODE_HALF_WIDTH = 1.2    # cm
NODE_HALF_HEIGHT = 0.9   # cm
MARGIN = 0.25            # extra margin around nodes
ROUTE_SPACING = 0.5      # spacing between parallel backward routes
ENTRY_OFFSET = 0.2       # vertical offset between backward edges entering one node


class TikZEdgeRouter:
    """Computes waypoints for edges; ``positions`` maps a node index to (x, y) and holds visible nodes only."""

    def __init__(self, positions, options):
        self.positions = positions
        self.options = options
        self.backward_edge_count = 0
        self.global_min_x = JAVA_DOUBLE_MAX_VALUE
        self.global_max_x = JAVA_DOUBLE_MIN_VALUE
        self.global_min_y = JAVA_DOUBLE_MAX_VALUE
        self.global_max_y = JAVA_DOUBLE_MIN_VALUE
        for pos in positions.values():
            self.global_min_x = jmin(self.global_min_x, pos[0] - NODE_HALF_WIDTH)
            self.global_max_x = jmax(self.global_max_x, pos[0] + NODE_HALF_WIDTH)
            self.global_min_y = jmin(self.global_min_y, pos[1] - NODE_HALF_HEIGHT)
            self.global_max_y = jmax(self.global_max_y, pos[1] + NODE_HALF_HEIGHT)

    @staticmethod
    def is_self_loop(frm, to):
        return frm == to

    isSelfLoop = is_self_loop

    def compute_waypoints(self, frm, to, all_nodes):
        """Waypoints [(x, y), ...] for the edge frm -> to; empty when a straight line is fine."""
        from_pos = self.positions.get(frm)
        to_pos = self.positions.get(to)
        if from_pos is None or to_pos is None:
            return []
        if frm == to:
            return []  # self loops use dedicated TikZ rendering
        if to_pos[0] < from_pos[0] - NODE_HALF_WIDTH:
            return self._route_backward_edge(from_pos, to_pos)
        if abs(to_pos[0] - from_pos[0]) < NODE_HALF_WIDTH:
            return self._route_vertical_edge(from_pos, to_pos)
        obstacles = self._find_obstacles(frm, to, from_pos, to_pos, all_nodes)
        if not obstacles:
            return []
        return self._route_forward_edge(from_pos, to_pos, obstacles)

    computeWaypoints = compute_waypoints

    @staticmethod
    def _route_vertical_edge(from_pos, to_pos):
        offset_x = NODE_HALF_WIDTH + MARGIN + 0.5
        return [(from_pos[0] + offset_x, from_pos[1]), (to_pos[0] + offset_x, to_pos[1])]

    def _route_backward_edge(self, from_pos, to_pos):
        # see _kb/12-interfaces-and-docs.md for the backward-edge routing (planarity) rationale
        self.backward_edge_count += 1
        count = self.backward_edge_count
        route_above = from_pos[1] >= to_pos[1]
        if route_above:
            route_y = self.global_max_y + MARGIN + (ROUTE_SPACING * count)
        else:
            route_y = self.global_min_y - MARGIN - (ROUTE_SPACING * count)
        exit_x = from_pos[0] + NODE_HALF_WIDTH + MARGIN + (count - 1) * 0.3
        entry_x = to_pos[0] - NODE_HALF_WIDTH - MARGIN - (count - 1) * 0.3
        target_y_offset = (count - 1) * ENTRY_OFFSET
        if not route_above:
            target_y_offset = -target_y_offset
        target_y = to_pos[1] + target_y_offset
        return [(exit_x, from_pos[1]), (exit_x, route_y), (entry_x, route_y), (entry_x, target_y)]

    def _route_forward_edge(self, from_pos, to_pos, obstacles):
        avg_y = 0.0
        for obs in obstacles:
            avg_y += self.positions[obs][1]
        avg_y /= len(obstacles)
        mid_y = (from_pos[1] + to_pos[1]) / 2.0
        route_above = mid_y >= avg_y
        route_y = avg_y + (NODE_HALF_HEIGHT + MARGIN + 0.15) * (1 if route_above else -1)

        min_x = JAVA_DOUBLE_MAX_VALUE
        max_x = JAVA_DOUBLE_MIN_VALUE
        for obs in obstacles:
            pos = self.positions[obs]
            min_x = jmin(min_x, pos[0] - NODE_HALF_WIDTH - MARGIN)
            max_x = jmax(max_x, pos[0] + NODE_HALF_WIDTH + MARGIN)

        entry_x = jmax(from_pos[0] + 0.3, min_x - 0.3)
        waypoints = [(entry_x, route_y)]
        exit_x = jmin(to_pos[0] - 0.3, max_x + 0.3)
        if exit_x > entry_x + 0.1:
            waypoints.append((exit_x, route_y))
        return waypoints

    def _find_obstacles(self, frm, to, from_pos, to_pos, all_nodes):
        obstacles = []
        for node in all_nodes:
            if node == frm or node == to:
                continue
            pos = self.positions.get(node)
            if pos is None:
                continue
            if _line_intersects_node(from_pos, to_pos, pos):
                obstacles.append(node)
        obstacles.sort(key=lambda node: self.positions[node][0])
        return obstacles

    @staticmethod
    def render_routed_edge(from_id, to_id, from_anchor, to_anchor, waypoints, prob, options):
        """TikZ code for an edge, straight or through orthogonal waypoints, with an optional probability label."""
        show = (options.show_routing_prob and not math.isnan(prob)
                and options.min_prob_to_show <= prob < 1.0 - options.min_prob_to_show)
        if not waypoints:
            if show:
                return '\\draw[conn] (%s%s) -- node[problabel] {%s} (%s%s);\n' % (
                    from_id, from_anchor, f2(prob), to_id, to_anchor)
            return '\\draw[conn] (%s%s) -- (%s%s);\n' % (from_id, from_anchor, to_id, to_anchor)
        parts = ['\\draw[conn] (%s%s)' % (from_id, from_anchor)]
        for wp in waypoints:
            parts.append(' -- (%s,%s)' % (f2(wp[0]), f2(wp[1])))
        if show:
            parts.append(' -- node[problabel] {%s} (%s%s);\n' % (f2(prob), to_id, to_anchor))
        else:
            parts.append(' -- (%s%s);\n' % (to_id, to_anchor))
        return ''.join(parts)

    renderRoutedEdge = render_routed_edge


def _line_intersects_node(line_start, line_end, center):
    half_w = NODE_HALF_WIDTH + MARGIN
    half_h = NODE_HALF_HEIGHT + MARGIN
    left = center[0] - half_w
    right = center[0] + half_w
    bottom = center[1] - half_h
    top = center[1] + half_h
    x1, y1 = line_start
    x2, y2 = line_end
    if (x1 < left and x2 < left) or (x1 > right and x2 > right):
        return False
    if (y1 < bottom and y2 < bottom) or (y1 > top and y2 > top):
        return False
    if _point_in_box(x1, y1, left, right, bottom, top) or _point_in_box(x2, y2, left, right, bottom, top):
        return True
    return (_segments_intersect(x1, y1, x2, y2, left, bottom, left, top)
            or _segments_intersect(x1, y1, x2, y2, right, bottom, right, top)
            or _segments_intersect(x1, y1, x2, y2, left, bottom, right, bottom)
            or _segments_intersect(x1, y1, x2, y2, left, top, right, top))


def _point_in_box(x, y, left, right, bottom, top):
    return left <= x <= right and bottom <= y <= top


def _segments_intersect(x1, y1, x2, y2, x3, y3, x4, y4):
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1e-10:
        return False
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom
    return 0 <= t <= 1 and 0 <= u <= 1
