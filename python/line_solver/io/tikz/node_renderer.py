"""
Renders queueing network nodes as TikZ code.

Port of jar/src/main/java/jline/io/tikz/TikZNodeRenderer.java: the preamble
with the node styles, and one renderer per node type.
"""

import math

from ._jcompat import escape_latex, f2, node_kind, sanitize_id

_STYLES = r"""\tikzset{
    % Queue: Rectangle buffer
    queue/.style={
        rectangle,
        draw=black,
        minimum width=1.1cm,
        minimum height=0.72cm,
        fill=white
    },
    % Server circle
    server/.style={
        circle,
        draw=black,
        minimum size=0.6cm,
        fill=white
    },
    % Delay: Vertical rectangle (infinite server)
    delay/.style={
        rectangle,
        draw=black,
        minimum width=0.8cm,
        minimum height=1.5cm,
        fill=green!10
    },
    % Source: Small circle (white)
    source/.style={
        circle,
        draw=black,
        minimum size=0.3cm,
        fill=white
    },
    % Sink: Small circle (black)
    sink/.style={
        circle,
        draw=black,
        minimum size=0.3cm,
        fill=black
    },
    % Fork: Diamond
    fork/.style={
        diamond,
        draw=black,
        minimum size=1cm,
        fill=orange!20,
        aspect=1.5
    },
    % Join: Diamond
    joinnode/.style={
        diamond,
        draw=black,
        minimum size=1cm,
        fill=purple!20,
        aspect=1.5
    },
    % Router: Hexagon
    router/.style={
        regular polygon,
        regular polygon sides=6,
        draw=black,
        minimum size=1cm,
        fill=cyan!10
    },
    % ClassSwitch: Trapezium
    classswitch/.style={
        trapezium,
        draw=black,
        trapezium left angle=70,
        trapezium right angle=110,
        minimum width=1.5cm,
        minimum height=0.8cm,
        fill=pink!20
    },
    % Cache: Stacked rectangle
    cache/.style={
        rectangle,
        draw=black,
        minimum width=1.5cm,
        minimum height=1.2cm,
        fill=gray!20
    },
    % Logger: Rectangle with lines
    logger/.style={
        rectangle,
        draw=black,
        minimum width=1.2cm,
        minimum height=0.8cm,
        fill=brown!10
    },
    % Place (Petri net): Circle
    place/.style={
        circle,
        draw=black,
        minimum size=0.8cm,
        fill=white
    },
    % Transition (Petri net): Rectangle
    transition/.style={
        rectangle,
        draw=black,
        minimum width=0.2cm,
        minimum height=1cm,
        fill=black
    },
    % Connection arrow
    conn/.style={
        ->,
        >=Stealth,
        thick
    },
    % Probability label
    problabel/.style={
        font=\footnotesize,
        midway,
        above,
        sloped
    },
    % Node name label
    nodelabel/.style={
        font=\small
    }
}
"""

# node kinds drawn as a single styled node plus the optional name label
_SIMPLE_STYLE = {'fork': 'fork', 'join': 'joinnode', 'router': 'router', 'classswitch': 'classswitch',
                 'logger': 'logger', 'place': 'place', 'transition': 'transition',
                 'generic': 'draw,circle,minimum size=0.8cm'}


def _java_int_str(value):
    """Render an int-valued option as Java's StringBuilder.append(int) would."""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


class TikZNodeRenderer:
    """Emits the preamble and the TikZ code for each node."""

    def __init__(self, options):
        self.options = options

    def get_preamble(self):
        return ('\\documentclass[tikz,border=' + _java_int_str(self.options.border_padding) + 'pt]{standalone}\n'
                '\\usepackage{tikz}\n'
                '\\usetikzlibrary{arrows.meta,positioning,shapes.geometric,shapes.misc,calc,decorations.pathreplacing}\n'
                '\n' + _STYLES)

    getPreamble = get_preamble

    def render_node(self, node, x, y):
        node_id = sanitize_id(node.get_name())
        kind = node_kind(node)
        if kind == 'queue':
            return self._render_queue(node, node_id, x, y)
        if kind == 'delay':
            code = '\\node[delay] (%s) at (%s,%s) {$\\infty$};\n' % (node_id, f2(x), f2(y))
            return code + self._name_label(node, node_id)
        if kind == 'source':
            code = '\\node[source] (%s) at (%s,%s) {};\n' % (node_id, f2(x), f2(y))
            code += self._name_label(node, node_id)
            return code + '\\draw[conn] ([xshift=-0.6cm]%s.west) -- (%s.west);\n' % (node_id, node_id)
        if kind == 'sink':
            code = '\\node[sink] (%s) at (%s,%s) {};\n' % (node_id, f2(x), f2(y))
            code += self._name_label(node, node_id)
            return code + '\\draw[conn] (%s.east) -- ([xshift=0.6cm]%s.east);\n' % (node_id, node_id)
        if kind == 'cache':
            code = '\\node[cache] (%s) at (%s,%s) {};\n' % (node_id, f2(x), f2(y))
            code += self._name_label(node, node_id)
            code += '\\draw ([yshift=-0.3cm]%s.north west) -- ([yshift=-0.3cm]%s.north east);\n' % (node_id, node_id)
            code += '\\draw (%s.west) -- (%s.east);\n' % (node_id, node_id)
            return code + '\\draw ([yshift=0.3cm]%s.south west) -- ([yshift=0.3cm]%s.south east);\n' % (node_id, node_id)
        code = '\\node[%s] (%s) at (%s,%s) {};\n' % (_SIMPLE_STYLE[kind], node_id, f2(x), f2(y))
        return code + self._name_label(node, node_id)

    renderNode = render_node

    def _name_label(self, node, node_id):
        if not self.options.show_node_names:
            return ''
        return '\\node[nodelabel,above=2pt of %s] {%s};\n' % (node_id, escape_latex(node.get_name()))

    def _render_queue(self, queue, node_id, x, y):
        code = '\\node[queue] (%s) at (%s,%s) {};\n' % (node_id, f2(x), f2(y))
        code += self._name_label(queue, node_id)
        if self.options.show_scheduling:
            sched = queue.get_sched_strategy()
            if sched is not None:
                name = getattr(sched, 'name', str(sched))
                # FCFSPRIO is MATLAB's alias of HOL (both are 11 there), so MATLAB labels it HOL
                code += '\\node[font=\\tiny,below=2pt of %s] {%s};\n' % (node_id, 'HOL' if name == 'FCFSPRIO' else name)
        if self.options.show_server_count:
            servers = queue.get_number_of_servers()
            if servers is None or math.isinf(float(servers)):
                label = '$\\infty$'
            else:
                label = str(int(servers))
            code += '\\node[server,anchor=west] (%s_server) at (%s.east) {%s};\n' % (node_id, node_id, label)
        return code

    def render_connection(self, from_node, to_node, prob):
        """Straight connection between two nodes, labelled with prob unless it is NaN or essentially 0 or 1."""
        from_id = sanitize_id(from_node.get_name())
        to_id = sanitize_id(to_node.get_name())
        fa = self._output_anchor(from_node)
        if self._show_prob(prob):
            return '\\draw[conn] (%s%s) -- node[problabel] {%s} (%s.west);\n' % (from_id, fa, f2(prob), to_id)
        return '\\draw[conn] (%s%s) -- (%s.west);\n' % (from_id, fa, to_id)

    renderConnection = render_connection

    def render_curved_connection(self, from_node, to_node, prob, bend_angle):
        """Curved connection for parallel paths, bending left for a positive angle."""
        from ._jcompat import jformat
        from_id = sanitize_id(from_node.get_name())
        to_id = sanitize_id(to_node.get_name())
        fa = self._output_anchor(from_node)
        bend = 'bend left' if bend_angle > 0 else 'bend right'
        angle = jformat(abs(bend_angle), 0)
        if self._show_prob(prob):
            return '\\draw[conn] (%s%s) to[%s=%s] node[problabel] {%s} (%s.west);\n' % (
                from_id, fa, bend, angle, f2(prob), to_id)
        return '\\draw[conn] (%s%s) to[%s=%s] (%s.west);\n' % (from_id, fa, bend, angle, to_id)

    renderCurvedConnection = render_curved_connection

    def _show_prob(self, prob):
        o = self.options
        return (o.show_routing_prob and not math.isnan(prob)
                and o.min_prob_to_show <= prob < 1.0 - o.min_prob_to_show)

    def _output_anchor(self, node):
        # a Delay draws no server circle; a Queue has one only when server counts are shown
        if node_kind(node) == 'queue' and self.options.show_server_count:
            return '_server.east'
        return '.east'
