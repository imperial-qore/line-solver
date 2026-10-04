"""
Exports a queueing network model as a TikZ diagram.

Port of jar/src/main/java/jline/io/tikz/TikZExporter.java: generates the
standalone LaTeX document, writes it to a .tex file, compiles it with
pdflatex, rasterizes it with pdftoppm and opens it in a PDF viewer.
"""

import os
import re
import shutil
import subprocess
import tempfile
from collections import deque

from ._jcompat import JAVA_DOUBLE_MAX_VALUE, JAVA_DOUBLE_MIN_VALUE, jformat, jmax, jmin, f2, node_kind, sanitize_id
from .edge_router import TikZEdgeRouter
from .layout import TikZLayoutEngine, _connmatrix
from .node_renderer import TikZNodeRenderer
from .options import TikZOptions
from .viewer import display_pdf


class TikZExporter:
    """Builds the TikZ document for a Network and drives the LaTeX toolchain."""

    def __init__(self, model, options=None):
        self.model = model
        self.options = TikZOptions.from_value(options)
        self.renderer = TikZNodeRenderer(self.options)
        self.layout_engine = TikZLayoutEngine(model, self.options)

    def get_options(self):
        return self.options

    getOptions = get_options

    # ------------------------------------------------------------------ text

    def generate_tikz(self):
        """Complete LaTeX document (standalone class) with the network drawn in TikZ."""
        return self._generate(environment=None)

    generateTikZ = generate_tikz

    def generate_environment_tikz(self, stage_names, transition_rates):
        """Network plus, to its right, the Markov chain of the environment stages with its transition rates."""
        return self._generate(environment=(list(stage_names), transition_rates))

    generateEnvironmentTikZ = generate_environment_tikz

    def _generate(self, environment):
        self.layout_engine.compute_layout()
        positions = self.layout_engine.get_all_positions()
        all_nodes = self.model.get_nodes()
        visible = self._filter_visible_nodes(all_nodes)
        collapsed = self._build_collapsed_connections(all_nodes, visible)

        out = [self.renderer.get_preamble(), '\n', '\\begin{document}\n']
        out.append('\\begin{tikzpicture}[scale=1.1]\n' if environment else '\\begin{tikzpicture}\n')
        out.append('\n')

        out.append('% Queueing Network\n' if environment else '% Nodes\n')
        for idx in visible:
            pos = positions.get(idx)
            if pos is not None:
                out.append(self.renderer.render_node(all_nodes[idx], pos[0], pos[1]))
        out.append('\n')

        out.append('% Network Connections\n' if environment else '% Connections\n')
        if _connmatrix(self.model) is not None:
            router = TikZEdgeRouter({idx: positions.get(idx) for idx in visible}, self.options)
            nan = float('nan')
            for frm in visible:
                for to in collapsed.get(frm, ()):
                    from_id = sanitize_id(all_nodes[frm].get_name())
                    to_id = sanitize_id(all_nodes[to].get_name())
                    if router.is_self_loop(frm, to):
                        out.append(self._render_self_loop(from_id, all_nodes[frm]))
                    else:
                        waypoints = router.compute_waypoints(frm, to, visible)
                        out.append(TikZEdgeRouter.render_routed_edge(
                            from_id, to_id, self._output_anchor(all_nodes[frm]), '.west', waypoints, nan, self.options))
        out.append('\n')

        if environment:
            out.append(self._render_environment_chain(positions, environment[0], environment[1]))
            out.append('\n')

        out.append('\\end{tikzpicture}\n')
        out.append('\\end{document}\n')
        return ''.join(out)

    @staticmethod
    def _render_environment_chain(positions, stage_names, rates):
        max_x = JAVA_DOUBLE_MIN_VALUE
        min_y = JAVA_DOUBLE_MAX_VALUE
        max_y = JAVA_DOUBLE_MIN_VALUE
        for pos in positions.values():
            max_x = jmax(max_x, pos[0])
            min_y = jmin(min_y, pos[1])
            max_y = jmax(max_y, pos[1])
        mc_x = max_x + 6.0
        n = len(stage_names)
        spacing = (max_y - min_y) / (n - 1) if n > 1 else 0.0
        if spacing < 1.5:
            spacing = 1.5
        center_y = (min_y + max_y) / 2.0

        out = ['% Environment Markov Chain\n']
        for i in range(n):
            y = center_y + (n - 1 - 2 * i) * spacing / 2.0
            out.append('\\node[draw,circle,minimum size=0.84cm,fill=blue!10] (mc_%d) at (%s,%s) {%s};\n'
                       % (i, f2(mc_x), f2(y), stage_names[i]))
        out.append('\n')
        out.append('% Environment Transitions\n')
        for i in range(n):
            for j in range(n):
                rate = float(rates[i][j])
                if i != j and rate > 0:
                    bend = 'bend left' if i < j else 'bend right'
                    out.append('\\draw[conn] (mc_%d) to[%s=50] node[midway,fill=white,font=\\footnotesize] {%s} (mc_%d);\n'
                               % (i, bend, jformat(rate, 1), j))
        return ''.join(out)

    # ------------------------------------------------------------ structure

    def _filter_visible_nodes(self, all_nodes):
        idx = list(range(len(all_nodes)))
        if not self.options.hide_auto_generated_nodes:
            return idx
        return [i for i in idx if not self._is_auto_generated(all_nodes[i])]

    @staticmethod
    def _is_auto_generated(node):
        if node_kind(node) != 'classswitch':
            return False
        name = node.get_name()
        return name.startswith('CS_') and '_to_' in name

    def _build_collapsed_connections(self, all_nodes, visible):
        """Visible successors of each visible node, walking through hidden nodes breadth first.

        Targets are in ascending node index, the connection-matrix row order, as the JAR (a node-index TreeSet)
        and the C++ port write them, and so MATLAB, which draws through the JAR.
        """
        result = {v: [] for v in visible}
        conn = _connmatrix(self.model)
        if conn is None:
            return result
        n = len(all_nodes)
        visible_set = set(visible)
        successors = {i: [j for j in range(n) if conn[i][j] > 0] for i in range(n)}
        for frm in visible:
            seen = set()
            queue = deque()
            for succ in successors[frm]:
                if succ not in seen:
                    queue.append(succ)
                    seen.add(succ)
            while queue:
                cur = queue.popleft()
                if cur in visible_set:
                    if cur not in result[frm]:
                        result[frm].append(cur)
                else:
                    for succ in successors[cur]:
                        if succ not in seen:
                            queue.append(succ)
                            seen.add(succ)
            result[frm].sort()
        return result

    def _render_self_loop(self, node_id, node):
        # see _kb/12-interfaces-and-docs.md for the self-loop rendering rationale
        # a Delay draws no server circle, so it takes the plain-node loop (as the JAR does since the same fix)
        if node_kind(node) == 'queue' and self.options.show_server_count:
            return ('\\draw[conn] (%s_server.north) -- ++(0,0.6) -- ++(-1.8,0) -- ++(0,-0.6) -- (%s.west);\n'
                    % (node_id, node_id))
        return '\\draw[conn] (%s.north) -- ++(0,0.5) -- ++(-0.8,0) -- ++(0,-0.5) -- (%s.west);\n' % (node_id, node_id)

    def _output_anchor(self, node):
        # a queue leaves from its server circle, which exists only when server counts are shown
        if node_kind(node) == 'queue' and self.options.show_server_count:
            return '_server.east'
        return '.east'

    # ---------------------------------------------------------------- files

    def export_to_file(self, file_path):
        """Write the TikZ document to file_path (UTF-8)."""
        with open(file_path, 'w', encoding='utf-8', newline='') as fh:
            fh.write(self.generate_tikz())

    exportToFile = export_to_file

    def export_to_pdf(self):
        """Compile the document with pdflatex in a fresh temporary directory and return the PDF path."""
        return self._compile(self.generate_tikz(), 'line-tikz-')

    exportToPDF = export_to_pdf

    def export_to_png(self, file_path, dpi=150):
        """Compile to PDF, then rasterize with pdftoppm into file_path at the given DPI."""
        _pdftoppm(self.export_to_pdf(), file_path, dpi)

    exportToPNG = export_to_png

    def export_environment_to_png(self, file_path, dpi, stage_names, transition_rates):
        """Compile the network plus environment chain diagram and rasterize it into file_path."""
        pdf = self._compile(self.generate_environment_tikz(stage_names, transition_rates), 'line-tikz-env-')
        _pdftoppm(pdf, file_path, dpi)

    exportEnvironmentToPNG = export_environment_to_png

    @staticmethod
    def _compile(tikz_code, prefix):
        if not is_pdflatex_available():
            raise RuntimeError('pdflatex is not available. Please install TeX Live or MiKTeX.')
        temp_dir = tempfile.mkdtemp(prefix=prefix)
        tex_file = os.path.join(temp_dir, 'network.tex')
        pdf_file = os.path.join(temp_dir, 'network.pdf')
        with open(tex_file, 'w', encoding='utf-8', newline='') as fh:
            fh.write(tikz_code)
        proc = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-output-directory=' + temp_dir, tex_file],
                              cwd=temp_dir, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
        if proc.returncode != 0 or not os.path.exists(pdf_file):
            raise RuntimeError('pdflatex compilation failed:\n' + proc.stdout.decode('utf-8', 'replace'))
        return pdf_file

    def display(self):
        """Compile to PDF and open it in a viewer; without pdflatex, save network-diagram.tex instead."""
        tex_path = 'network-diagram.tex'
        if not is_pdflatex_available():
            try:
                self.export_to_file(tex_path)
                print('pdflatex not found. TikZ code saved to: ' + tex_path)
                print('Compile manually with: pdflatex ' + tex_path)
            except OSError as exc:
                print('Failed to save TikZ file: ' + str(exc))
            return None
        try:
            pdf = self.export_to_pdf()
            display_pdf(pdf, self.model.get_name())
            return pdf
        except (OSError, RuntimeError) as exc:
            print('Failed to generate visualization: ' + str(exc))
            try:
                self.export_to_file(tex_path)
                print('TikZ code saved to: ' + tex_path)
            except OSError as exc2:
                print('Failed to save TikZ file: ' + str(exc2))
            return None


def is_pdflatex_available():
    """True when a working pdflatex is on PATH."""
    if shutil.which('pdflatex') is None:
        return False
    try:
        return subprocess.run(['pdflatex', '--version'], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL).returncode == 0
    except OSError:
        return False


def is_pdftoppm_available():
    """True when pdftoppm (poppler-utils) is on PATH."""
    return shutil.which('pdftoppm') is not None


def _pdftoppm(pdf_file, png_path, dpi):
    stem = re.sub(r'\.png$', '', png_path)
    try:
        proc = subprocess.run(['pdftoppm', '-png', '-r', str(int(dpi)), '-singlefile', os.path.abspath(pdf_file), stem],
                              stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except FileNotFoundError as exc:
        raise RuntimeError('pdftoppm is not available. Please install poppler-utils.') from exc
    if proc.returncode != 0:
        raise RuntimeError('pdftoppm conversion failed:\n' + proc.stdout.decode('utf-8', 'replace'))
