#!/usr/bin/env python3
"""Derive the Jupyter example notebooks from the ``.py`` example scripts.

The notebooks under ``python/examples/`` used to be hand-authored beside the
scripts, so the two drifted: model parameters diverged, notebooks called
solvers the scripts never called, cells lost their newlines and silently
collapsed into a single comment, and 206 scripts had no notebook at all.  This
generator makes the notebook a DERIVED artifact: every code cell it emits comes
verbatim from the script, so the only way a notebook can disagree with its
script is for someone to hand-edit the notebook and skip the regeneration.

Cell boundaries come from ``# %%`` / ``# %% [markdown]`` markers when the script
carries them (jupytext percent format), and otherwise from the script's own
top-level structure: the module docstring becomes the title, imports group into
a setup cell, each ``def``/``class`` gets its own cell, and the remaining
statements split on blank lines exactly where the author put them.

An ``if __name__ == '__main__':`` guard is unwrapped so the notebook actually
runs the demo, and ``sys.exit(main())`` becomes ``main()`` so the kernel is not
killed.  Both rewrites are checked: :func:`verify` reparses the emitted cells
and compares the statement-level AST against the script's, so a split that
would change what runs is a hard error rather than a silent divergence.

Usage::

    python3 tools/gen_notebooks.py              # regenerate every notebook
    python3 tools/gen_notebooks.py --check      # fail if any is out of date
    python3 tools/gen_notebooks.py examples/basic/openQN/oqn_oneline.py
"""

import argparse
import ast
import json
import os
import re
import sys

EXAMPLES_DIR = 'examples'

# The root quickstart lives outside examples/ but carries a notebook twin that
# test_quickstart_scripts.py exercises, so it is generated with the rest.
EXTRA_EXAMPLES = ('mm1.py',)

# Directories and files that are not runnable examples and get no notebook.
# test_gallery/ holds pytest cases, not examples; the remaining skips are the
# helper modules the generator's own runnable-body test filters out.
SKIP_DIRS = ('examples/test_gallery',)
SKIP_NAMES = ('__init__.py',)

PERCENT_RE = re.compile(r'^# %%(?:\s+\[(?P<kind>\w+)\])?\s*$')
SHEBANG_RE = re.compile(r'^#!')
CODING_RE = re.compile(r'^#.*coding[:=]')

KERNELSPEC = {
    'display_name': 'Python 3',
    'language': 'python',
    'name': 'python3',
}

LANGUAGE_INFO = {
    'codemirror_mode': {'name': 'ipython', 'version': 3},
    'file_extension': '.py',
    'mimetype': 'text/x-python',
    'name': 'python',
    'nbconvert_exporter': 'python',
    'pygments_lexer': 'ipython3',
}

# Bound to the notebook's own directory so that a script locating a data file
# through __file__ resolves it exactly as it does when run as a script.
FILE_SHIM = """\
# A Jupyter kernel defines no __file__; bind it to this notebook's directory so
# that path lookups in the example resolve as they do for the .py script.
import os
__file__ = os.path.join(os.getcwd(), {stem!r})
"""


class GenError(Exception):
    """A script the generator refuses to convert rather than convert wrongly."""


# --------------------------------------------------------------------------
# discovery
# --------------------------------------------------------------------------

def is_example(path):
    """True when ``path`` is a runnable example that should carry a notebook."""
    posix = path.replace(os.sep, '/')
    if any(posix.startswith(d + '/') for d in SKIP_DIRS):
        return False
    if os.path.basename(path) in SKIP_NAMES:
        return False
    try:
        tree = ast.parse(open(path).read())
    except SyntaxError:
        return False
    return bool(runnable_body(tree))


def runnable_body(tree):
    """Top-level statements that do something when the module runs as a script.

    A file holding only imports, defs and classes is a helper module, not an
    example, and gets no notebook.
    """
    out = []
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom,
                             ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            continue
        if is_docstring(node):
            continue
        out.append(node)
    return out


def find_examples(root=EXAMPLES_DIR):
    paths = [p for p in EXTRA_EXAMPLES if os.path.exists(p) and is_example(p)]
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for name in sorted(filenames):
            if not name.endswith('.py'):
                continue
            path = os.path.join(dirpath, name)
            if is_example(path):
                paths.append(path)
    return paths


# --------------------------------------------------------------------------
# source helpers
# --------------------------------------------------------------------------

def is_docstring(node):
    return (isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str))


def is_main_guard(node):
    """True for ``if __name__ == '__main__':`` in any of its spellings."""
    if not isinstance(node, ast.If):
        return False
    test = node.test
    if not isinstance(test, ast.Compare) or len(test.comparators) != 1:
        return False
    if not isinstance(test.ops[0], ast.Eq):
        return False
    left, right = test.left, test.comparators[0]
    named = (isinstance(left, ast.Name) and left.id == '__name__')
    valued = (isinstance(right, ast.Constant) and right.value == '__main__')
    return named and valued


def strip_header(lines):
    """Drop a shebang and an encoding declaration; they mean nothing in a cell."""
    i = 0
    while i < len(lines) and (SHEBANG_RE.match(lines[i]) or CODING_RE.match(lines[i])):
        i += 1
    return lines[i:], i


def dedent_block(lines, width):
    """Remove exactly ``width`` columns of indentation from every non-blank line.

    Unlike ``textwrap.dedent`` this never touches a line that is shorter than
    the indent (a continuation line inside a triple-quoted string), so a block
    it cannot dedent safely is caught by the caller's AST comparison instead of
    being silently rewritten.
    """
    out = []
    pad = ' ' * width
    for line in lines:
        if not line.strip():
            out.append('\n' if line.endswith('\n') else '')
        elif line.startswith(pad):
            out.append(line[width:])
        elif line.startswith('\t'):
            raise GenError('tab-indented block cannot be dedented reliably')
        else:
            out.append(line)
    return out


def comment_prefix(lines, start, floor):
    """Index of the first line of the comment block attached above ``start``.

    Comments and blank lines between the previous statement and this one belong
    to this statement, which is what keeps a section header comment glued to the
    code it introduces.
    """
    i = start
    while i > floor:
        stripped = lines[i - 1].strip()
        if stripped.startswith('#') or not stripped:
            i -= 1
        else:
            break
    # blank lines immediately above the statement are separators, not comment
    while i < start and not lines[i].strip():
        i += 1
    return i


def statement_chunks(stmts, lines, floor):
    """Split ``stmts`` into ``(text, node, blank_before)`` source chunks.

    Statements sharing a physical line (``a(); b()``) stay in one chunk, since
    the chunk is a slice of source lines and emitting one per statement would
    repeat the whole line once per semicolon.
    """
    chunks = []
    heads = []
    prev_end = floor
    for node in stmts:
        end = node.end_lineno
        for deco in getattr(node, 'decorator_list', []):
            end = max(end, deco.end_lineno)
        start = min([node.lineno] + [d.lineno for d in
                                     getattr(node, 'decorator_list', [])]) - 1
        if chunks and start < prev_end:
            # continues a line already covered by the previous chunk
            prev_end = max(prev_end, end)
            chunks[-1] = (''.join(lines[heads[-1]:prev_end]).rstrip('\n'),
                          chunks[-1][1], chunks[-1][2])
            continue
        head = comment_prefix(lines, start, prev_end)
        blank_before = any(not lines[i].strip() for i in range(prev_end, head))
        chunks.append((''.join(lines[head:end]).rstrip('\n'), node, blank_before))
        heads.append(head)
        prev_end = end
    return chunks


# --------------------------------------------------------------------------
# splitting
# --------------------------------------------------------------------------

def title_from_stem(stem):
    """Fallback title: the file's own name, not prose invented around it."""
    return '`%s`' % stem


def docstring_markdown(doc, stem):
    """Turn a module docstring into the notebook's title cell."""
    body = [ln.rstrip() for ln in (doc or '').strip().splitlines()]
    while body and not body[0]:
        body.pop(0)
    if not body:
        return '# ' + title_from_stem(stem)
    title = body[0].rstrip('.')
    rest = [ln for ln in body[1:]]
    while rest and not rest[0]:
        rest.pop(0)
    out = '# ' + title
    if rest:
        out += '\n\n' + '\n'.join(rest).rstrip()
    return out


def split_percent(src, stem):
    """Split a jupytext percent-format script on its ``# %%`` markers."""
    lines = src.splitlines(keepends=True)
    lines, _ = strip_header(lines)
    cells = []
    kind = 'code'
    buf = []

    def flush():
        text = ''.join(buf).strip('\n')
        if not text:
            return
        if kind == 'markdown':
            stripped = []
            for ln in text.splitlines():
                stripped.append(ln[2:] if ln.startswith('# ') else
                                ln[1:] if ln.startswith('#') else ln)
            cells.append(('markdown', '\n'.join(stripped).strip('\n')))
        else:
            cells.append(('code', text))

    seen_marker = False
    for line in lines:
        m = PERCENT_RE.match(line.rstrip('\n'))
        if m:
            if seen_marker:
                flush()
            buf = []
            kind = 'markdown' if (m.group('kind') or '').lower() == 'markdown' else 'code'
            seen_marker = True
            continue
        buf.append(line)
    flush()
    return cells


def split_auto(src, stem):
    """Split a plain script on its own top-level structure."""
    lines = src.splitlines(keepends=True)
    lines, _ = strip_header(lines)
    tree = ast.parse(''.join(lines))
    cells = []

    body = list(tree.body)
    doc = None
    if body and is_docstring(body[0]):
        doc = body[0].value.value
        body = body[1:]
    cells.append(('markdown', docstring_markdown(doc, stem)))

    floor = 0
    if doc is not None:
        floor = tree.body[0].end_lineno

    groups = []          # list of (stmts, dedented_lines, floor)
    plain = []
    for node in body:
        if is_main_guard(node):
            if plain:
                groups.append((plain, lines, floor))
                floor = plain[-1].end_lineno
                plain = []
            groups.append(unwrap_main(node, lines))
            floor = node.end_lineno
        else:
            plain.append(node)
    if plain:
        groups.append((plain, lines, floor))

    for stmts, glines, gfloor in groups:
        cells.extend(group_cells(stmts, glines, gfloor))
    return cells


def unwrap_main(node, lines):
    """Return the guard's body as top-level statements over dedented lines."""
    inner = list(node.body)
    if not inner:
        raise GenError('empty __main__ guard')
    start = inner[0].lineno - 1
    end = inner[-1].end_lineno
    width = inner[0].col_offset
    block = dedent_block(lines[start:end], width)
    text = ''.join(block)
    try:
        sub = ast.parse(text)
    except SyntaxError as exc:
        raise GenError('cannot dedent __main__ guard: %s' % exc)
    if [ast.dump(s) for s in sub.body] != [ast.dump(s) for s in inner]:
        raise GenError('dedenting the __main__ guard changed its meaning')
    # `sub` is parsed from `block` alone, so its line numbers index `block`.
    return sub.body, block, 0


def group_cells(stmts, lines, floor):
    """Group statements into cells: imports together, defs alone, else by blanks."""
    chunks = statement_chunks(stmts, lines, floor)
    cells = []
    buf = []
    prev = None
    for text, node, blank_before in chunks:
        is_def = isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        is_imp = isinstance(node, (ast.Import, ast.ImportFrom))
        prev_def = isinstance(prev, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        prev_imp = isinstance(prev, (ast.Import, ast.ImportFrom))
        if prev is None:
            new = False
        elif is_def or prev_def:
            new = True
        elif is_imp != prev_imp:
            new = True
        else:
            new = blank_before
        if new and buf:
            cells.append(('code', '\n'.join(buf)))
            buf = []
        buf.append(text)
        prev = node
    if buf:
        cells.append(('code', '\n'.join(buf)))
    return cells


def rewrite_exit(cells):
    """Turn ``sys.exit(main())`` into ``main()`` so the kernel survives the cell.

    Only the unwrapped guard's own exit is rewritten: it exists to give the
    shell a status code, which a notebook has no use for, and SystemExit would
    otherwise be reported as a failed cell.
    """
    out = []
    for kind, text in cells:
        if kind == 'code':
            text = re.sub(r'^sys\.exit\((.+)\)\s*$', r'\1', text, flags=re.M)
            text = re.sub(r'^sys\.exit\(\s*\)\s*$', 'pass', text, flags=re.M)
        out.append((kind, text))
    return out


# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------

def canonical_script_ast(src):
    """Statement dumps of the script as the notebook is expected to run it."""
    tree = ast.parse(src)
    body = list(tree.body)
    if body and is_docstring(body[0]):
        body = body[1:]
    flat = []
    for node in body:
        if is_main_guard(node):
            flat.extend(node.body)
        else:
            flat.append(node)
    return [ast.dump(strip_exit(node)) for node in flat]


def strip_exit(node):
    """Apply the same exit rewrites :func:`rewrite_exit` makes to the cells."""
    if (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
            and dotted(node.value.func) == 'sys.exit'
            and not node.value.keywords):
        if len(node.value.args) == 1:
            return ast.Expr(value=node.value.args[0])
        if not node.value.args:
            return ast.Pass()
    return node


def dotted(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        base = dotted(node.value)
        return base + '.' + node.attr if base else None
    return None


def verify(src, cells):
    """Assert the notebook's code is the script's code, statement for statement."""
    code = '\n'.join(text for kind, text in cells if kind == 'code')
    try:
        tree = ast.parse(code)
    except SyntaxError as exc:
        raise GenError('generated cells do not parse: %s' % exc)
    body = list(tree.body)
    # the module docstring rides in the title cell, so either side may lack it
    if body and is_docstring(body[0]):
        body = body[1:]
    got = [ast.dump(n) for n in body]
    want = canonical_script_ast(src)
    if got != want:
        for i, (a, b) in enumerate(zip(got, want)):
            if a != b:
                raise GenError('statement %d differs:\n  notebook: %s\n  script:   %s'
                               % (i, a[:200], b[:200]))
        raise GenError('notebook has %d top-level statements, script has %d'
                       % (len(got), len(want)))


# --------------------------------------------------------------------------
# notebook assembly
# --------------------------------------------------------------------------

def build_cells(path):
    src = open(path).read()
    stem = os.path.splitext(os.path.basename(path))[0]
    if any(PERCENT_RE.match(ln) for ln in src.splitlines()):
        cells = split_percent(src, stem)
        title = docstring_markdown(ast.get_docstring(ast.parse(src)), stem)
        cells.insert(0, ('markdown', title))
    else:
        cells = split_auto(src, stem)
    cells = rewrite_exit(cells)
    verify(src, cells)
    if re.search(r'\b__file__\b', src):
        cells.insert(1, ('code', FILE_SHIM.format(stem=stem + '.py').rstrip('\n')))
    return cells


def to_notebook(cells):
    out = []
    for i, (kind, text) in enumerate(cells):
        cell = {
            'cell_type': kind,
            'id': 'cell-%d' % i,
            'metadata': {},
            'source': as_source(text),
        }
        if kind == 'code':
            cell['execution_count'] = None
            cell['outputs'] = []
        out.append(cell)
    return {
        'cells': out,
        'metadata': {'kernelspec': KERNELSPEC, 'language_info': LANGUAGE_INFO},
        'nbformat': 4,
        'nbformat_minor': 5,
    }


def as_source(text):
    """nbformat line list: every line keeps its newline except the last."""
    lines = text.split('\n')
    return [ln + '\n' for ln in lines[:-1]] + [lines[-1]]


def render(path):
    return json.dumps(to_notebook(build_cells(path)),
                      indent=1, ensure_ascii=False) + '\n'


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('paths', nargs='*',
                    help='example .py files (default: every example)')
    ap.add_argument('--check', action='store_true',
                    help='report notebooks that are out of date, write nothing')
    ap.add_argument('-q', '--quiet', action='store_true')
    args = ap.parse_args(argv)

    paths = args.paths or find_examples()
    stale, failed, wrote = [], [], 0
    for path in paths:
        nbpath = os.path.splitext(path)[0] + '.ipynb'
        try:
            text = render(path)
        except (GenError, SyntaxError) as exc:
            failed.append((path, exc))
            continue
        old = open(nbpath).read() if os.path.exists(nbpath) else None
        if old == text:
            continue
        if args.check:
            stale.append(nbpath)
        else:
            with open(nbpath, 'w') as fh:
                fh.write(text)
            wrote += 1

    for path, exc in failed:
        print('FAILED %s: %s' % (path, exc), file=sys.stderr)
    if args.check:
        for nbpath in stale:
            print('OUT OF DATE %s' % nbpath, file=sys.stderr)
        if stale or failed:
            print('%d notebook(s) out of date, %d script(s) failed; run '
                  'python3 tools/gen_notebooks.py' % (len(stale), len(failed)),
                  file=sys.stderr)
            return 1
        if not args.quiet:
            print('%d notebook(s) up to date' % len(paths))
        return 0
    if not args.quiet:
        print('%d example(s), %d notebook(s) written, %d failed'
              % (len(paths), wrote, len(failed)))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
