"""Shared harness for the notebook example suites.

`test_all_basic_examples.py`, `test_all_advanced_examples.py`,
`test_all_gettingstarted.py`, `test_all_gallery_examples.py` and
`test_all_toplevel_examples.py` all execute `python/examples/**/*.ipynb` and
compare the numbers each cell printed against a stored baseline.  The harness
used to be copy-pasted into every one of them, so the output filter -- which is
the comparison's real contract, see `_kb/08-build-and-test.md` -- had to be
edited in three places at once or the suites silently disagreed about what
counts as a result.  It lives here once instead.

A suite is now a subclass that names the directories it owns:

    class TestBasicExamples(NotebookSuite):
        roots = ['basic']

    attach(TestBasicExamples)

`roots` are collected RECURSIVELY, so a new example subfolder is picked up by
the suite that owns its parent without anyone editing a list.  The hardcoded
`enabled_subfolders` this replaces is what left 72 notebooks -- the whole of
`gallery/`, `basic/workflowModels`, `advanced/agentModels`, `advanced/fcRegion`
and `gettingstarted` tut11-tut13 -- in no suite at all, which is how 15
unparseable notebooks survived in the tree.
"""

import fcntl
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import traceback
import unittest
import warnings

import nbformat
import numpy as np
from nbconvert.preprocessors import ExecutePreprocessor

_REPO_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, _REPO_PYTHON)

# Each notebook runs with the kernel's cwd set to its OWN directory, so the
# in-tree package is not on the kernel's path and `import line_solver` binds
# whatever copy site-packages holds -- the retired line-solver-wrapper installs
# one under that very name, and cells using anything it predates die with a
# NameError that reads like a missing export. The kernel inherits this process's
# environment, and PYTHONPATH sits ahead of site-packages, so it pins the
# in-tree package for every notebook regardless of cwd.
os.environ['PYTHONPATH'] = os.pathsep.join(
    [_REPO_PYTHON] + [p for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p])

# Install a kernel for the current Python interpreter to ensure notebooks run
# with the same environment. The name carries a digest of that interpreter's
# path: the spec lives in the user-global jupyter directory, so a single fixed
# name is CLOBBERED by any concurrent run under a different interpreter (a
# system python3 run while the poetry venv suite is mid-flight repoints the
# venv suite's kernel at the system one, silently, mid-run).
KERNEL_NAME = 'line_solver_test_kernel_' + hashlib.sha1(
    sys.executable.encode()).hexdigest()[:8]

# THE INSTALL IS SERIALIZED ACROSS PROCESSES, because under `pytest -n` every
# xdist worker imports this module at once and they all install the SAME spec
# name. `install_kernel_spec` rmtree's the destination before it copies into it,
# so two concurrent installs open a window where the spec directory does not
# exist -- a third worker resolving the kernel then dies with NoSuchKernel, or
# the install itself fails and that worker silently falls back to the `python3`
# kernel, which is the system interpreter and not the venv the suite is pinned
# to. The lock is per-interpreter (the same digest that names the kernel), taken
# for the install only, and the spec is idempotent, so the first worker through
# does the work and the rest re-do it harmlessly, one at a time.
_KERNEL_LOCK = os.path.join(tempfile.gettempdir(), KERNEL_NAME + '.lock')
try:
    with open(_KERNEL_LOCK, 'w') as _lockf:
        fcntl.flock(_lockf, fcntl.LOCK_EX)
        subprocess.run(
            [sys.executable, '-m', 'ipykernel', 'install', '--user',
             '--name', KERNEL_NAME, '--display-name', 'LINE Test Kernel'],
            check=True, capture_output=True
        )
except (subprocess.CalledProcessError, OSError):
    # Fall back to default python3 kernel if installation fails
    KERNEL_NAME = 'python3'

# THE SAME CONTENTION ALLOWANCE THE PARITY ROWS GET, for the same reason and
# from the same measurement: under `pytest -n auto` a notebook shares its host
# with fifteen siblings, each running a kernel of its own, so its wall clock is
# several times its measured idle cost while nothing about it is unhealthy.
# `advanced/largeScale/largescale_tbi` was killed on the first parallel run
# (picard01, 16 workers, 2026-08-28) against the 1200 s default. Capped at 4x so
# a notebook that genuinely hangs is still cut off inside the phase. See the
# twin in tests/parity/test_parity_static.py.
_XDIST_WORKERS = int(os.environ.get('PYTEST_XDIST_WORKER_COUNT') or 1)
CONTENTION = max(1, min(_XDIST_WORKERS, 4))

# THE CELL BUDGET IS NOT THE ONLY DEADLINE, and the other one had no allowance.
# `wait_for_ready` gives the kernel `startup_timeout` to answer kernel_info
# before nbclient raises `RuntimeError: Kernel didn't respond in N seconds`, and
# that clock runs BEFORE the first cell, so `timeout` above never covers it.
# nbclient defaults it to 60 s. Measured idle here (picard09, otherwise idle):
# 2.1-3.6 s, so the default is already ~20x the cost -- and phase 5 still lost
# `basic/mixedQN/mqn_singleserver_fcfs` to it on picard05 on 2026-09-09, where 8
# workers each held a kernel plus a JMT java and an LDES subprocess on 16 cores.
# The notebook never ran a line: the failure was the boot, not the model.
# 120 s carries the same CONTENTION factor the cell budget gets, so `-n auto`
# allows up to 480 s while a serial run keeps a tight 120 s. It stays a real
# ceiling -- a kernel that DIED is reported at once by the liveness check, not
# at the deadline.
KERNEL_STARTUP = 120

TOLERANCE = 1e-6
RELATIVE_TOLERANCE = 1e-5


class NotebookSuite(unittest.TestCase):
    """Executes every notebook under `roots` and diffs its numbers vs baseline."""

    # pytest collects every unittest.TestCase it imports, name notwithstanding;
    # each concrete suite sets this back to True.
    __test__ = False

    notebooks_dir = 'examples'
    baselines_dir = 'tests/regression'

    #: directory prefixes under `notebooks_dir`, collected recursively.
    #: '' means the notebooks_dir itself, non-recursively.
    roots = []

    #: seconds a single notebook may run before the kernel is killed
    timeout = 1200

    #: Per-notebook overrides of `timeout`, keyed by the path under
    #: `notebooks_dir` -- the same string `run_notebook` is handed.
    #:
    #: A NOTEBOOK BELONGS HERE ONCE ITS IDLE COST IS A LARGE FRACTION OF THE
    #: DEFAULT, because past that point the budget stops being a hang detector
    #: and becomes a contention detector: the suite shares its host with
    #: whatever else the cluster is running, so the kernel is killed mid-cell
    #: and the run reports a FAILURE against code that is fine. This is the
    #: notebook twin of `SLOW_ROWS` in tests/parity/test_parity_static.py,
    #: widened on 2026-08-27 for exactly this reason.
    SLOW_NOTEBOOKS = {
        # SSA at 100k samples over Cache(1000 items, capacity 50, LRU) read
        # through Zipf(1.4, 1000). Measured ALONE at 840 s on 2026-08-28
        # (picard09, otherwise idle) -- 70% of the 1200 s default, and its RSS
        # climbs past 3 GB while it runs. It died on picard03 on 2026-08-27,
        # where phase 5 shared the node with the SolverNN datagen cluster; the
        # parity twin of this same example, `tut06_cache_lru_zipf[PYTHON]`,
        # failed in that same run for the same reason and was widened in
        # 62177fd1e. 3600 s is ~4x the idle cost, against the ~3x that twin got.
        'gettingstarted/tut06_cache_lru_zipf.ipynb': 3600,
    }

    def setUp(self):
        """Set up test environment."""
        self.working_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        self.baseline_warnings = []

    def run_notebook(self, filename):
        os.chdir(self.working_dir)
        filepath = os.path.join(self.notebooks_dir, filename)
        print(f'\n{"="*60}', flush=True)
        print(f'RUNNING NOTEBOOK: {filename}', flush=True)
        print(f'Working directory: {self.working_dir}', flush=True)
        print(f'Full path: {filepath}', flush=True)
        print(f'{"="*60}', flush=True)

        with open(filepath) as f:
            nb = nbformat.read(f, as_version=4)
            budget = self.SLOW_NOTEBOOKS.get(filename, self.timeout)
            # `timeout = None` means a suite has deliberately opted out of any
            # deadline (the basic models are cheap; a hang there is a defect and
            # should hang visibly), and None has no contention allowance to give.
            if budget is not None:
                budget *= CONTENTION
            ep = ExecutePreprocessor(timeout=budget, kernel_name=KERNEL_NAME,
                                     startup_timeout=KERNEL_STARTUP * CONTENTION)
            try:
                ep.preprocess(nb, {'metadata': {'path': os.path.dirname(filepath)}})
                print(f'✓ Successfully executed notebook: {filename}', flush=True)
            except Exception as e:
                print(f"\n{'='*60}", flush=True)
                print(f"NOTEBOOK EXECUTION FAILED: {filename}", flush=True)
                print(f"{'='*60}", flush=True)
                print(f"Error Type: {type(e).__name__}", flush=True)
                print(f"Error Message: {str(e)}", flush=True)
                print(f"\nFull Traceback:", flush=True)
                print(traceback.format_exc(), flush=True)
                print(f"{'='*60}", flush=True)
                raise

        for cell_idx, cell in enumerate(nb.cells):
            if cell.cell_type != 'code':
                continue
            for output in cell.get('outputs', []):
                if output.output_type != 'error':
                    continue
                error_name = output.get("ename", "")
                error_value = output.get("evalue", "")
                traceback_lines = output.get("traceback", [])
                cell_source = cell.get('source', '')

                print(f"\n{'='*60}", flush=True)
                print(f"CELL EXECUTION ERROR: {filename}, Cell {cell_idx}", flush=True)
                print(f"{'='*60}", flush=True)
                print(f"Error Name: {error_name}", flush=True)
                print(f"Error Value: {error_value}", flush=True)
                print(f"\nCell Source Code:", flush=True)
                print(f"{'-'*40}", flush=True)
                print(f"{cell_source}", flush=True)
                print(f"{'-'*40}", flush=True)
                print(f"\nFull Traceback:", flush=True)
                for line in traceback_lines:
                    print(line, flush=True)
                print(f"{'='*60}", flush=True)

                self.fail(
                    f"Notebook {filename} failed.\n"
                    f"Error Name: {error_name}\n"
                    f"Error Value: {error_value}\n"
                    f"See detailed output above for full context."
                )

        self._compare_with_baseline(filename, nb)

        return nb

    def _compare_with_baseline(self, filename, nb):
        """Compare notebook results with saved regression data."""
        notebook_dir = os.path.dirname(filename)
        notebook_name = os.path.basename(filename).replace('.ipynb', '_regression.json')
        regression_path = os.path.join(self.working_dir, self.baselines_dir,
                                       self.notebooks_dir, notebook_dir, notebook_name)

        if not os.path.exists(regression_path):
            warning_msg = (f"No regression data found for {filename}. "
                           f"Run gen_all_regression.py to generate regression data.")
            warnings.warn(warning_msg, UserWarning)
            self.baseline_warnings.append(warning_msg)
            return

        try:
            with open(regression_path, 'r') as f:
                regression_results = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            warning_msg = f"Failed to load regression data for {filename}: {e}"
            warnings.warn(warning_msg, UserWarning)
            self.baseline_warnings.append(warning_msg)
            return

        regression_results = self._flatten_baseline(regression_results)
        current_results = self._extract_numerical_outputs(nb)

        differences = self._compare_results(regression_results, current_results, filename)

        if differences:
            warning_msg = (f"Numerical differences detected in {filename}:\n"
                           + "\n".join(differences))
            warnings.warn(warning_msg, UserWarning)
            self.baseline_warnings.append(warning_msg)

    def _flatten_baseline(self, regression_results):
        """
        Reduce a stored baseline to the {cell_i: [numbers]} map used for comparison.

        gen_all_regression.py writes the comprehensive v2.0 schema
        ({'cell_content', 'numerical_summary', 'metadata'}); older baselines are
        already the flat map. The per-cell numbers are re-extracted here from the
        raw output text held in 'cell_content' with this class's own extractor, so
        generator and comparator cannot drift apart in what they count as a value
        (in particular 'numerical_summary' also folds in literals from the cell
        source, which the notebook-side extractor never sees).
        """
        if 'cell_content' not in regression_results:
            return regression_results

        flattened = {}
        for cell_key, cell_data in regression_results['cell_content'].items():
            if cell_data.get('cell_type') != 'code':
                continue
            cell_results = []
            for output in cell_data.get('outputs', []):
                content = output.get('content', {})
                output_type = output.get('output_type')
                if output_type in ('execute_result', 'display_data'):
                    text_output = content.get('text/plain')
                elif output_type == 'stream' and content.get('stream_name') == 'stdout':
                    text_output = content.get('text')
                else:
                    continue
                cell_results.extend(self._extract_numbers_from_text(text_output))
            if cell_results:
                flattened[cell_key] = cell_results
        return flattened

    def _extract_numerical_outputs(self, nb):
        """Extract numerical values from notebook outputs."""
        results = {}

        for cell_idx, cell in enumerate(nb.cells):
            if cell.cell_type == 'code':
                cell_results = []

                for output in cell.get('outputs', []):
                    if output.output_type in ['execute_result', 'display_data']:
                        if 'data' in output and 'text/plain' in output['data']:
                            text_output = output['data']['text/plain']
                            numerical_values = self._extract_numbers_from_text(text_output)
                            if numerical_values:
                                cell_results.extend(numerical_values)

                    elif output.output_type == 'stream' and output.name == 'stdout':
                        text_output = output['text']
                        numerical_values = self._extract_numbers_from_text(text_output)
                        if numerical_values:
                            cell_results.extend(numerical_values)

                if cell_results:
                    results[f'cell_{cell_idx}'] = cell_results

        return results

    def _extract_numbers_from_text(self, text):
        """Extract numerical values from text output, filtering runtimes and timestamps."""
        if not isinstance(text, str):
            return []

        text = self._filter_runtime_and_timestamps(text)

        number_pattern = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'

        matches = re.findall(number_pattern, text)

        numbers = []
        for match in matches:
            try:
                num = float(match)
                if not (np.isnan(num) or np.isinf(num)):
                    numbers.append(num)
            except (ValueError, OverflowError):
                continue

        return numbers

    def _filter_runtime_and_timestamps(self, text):
        """Filter out runtime and timestamp information from text."""
        if not isinstance(text, str):
            return text

        patterns_to_remove = [
            r'completed in \d+\.?\d*s\.?',
            r'\d{4}-\d{2}-\d{2}[\sT]\d{2}:\d{2}:\d{2}',
            r'\d{2}:\d{2}:\d{2}',
            r'\d{2}/\d{2}/\d{4}',
            r'\d{4}/\d{2}/\d{2}',
            r'/tmp/[^/\s]*\d{10,}[^/\s]*',
            r'/workspace/[^/\s]*\d{10,}[^/\s]*',
            r'JMT Model: [^\n]*',
            # Solver banner: "<NAME> analysis [method: ...; type: ...;
            # lang: ...; env: ...] completed in <t>s." with its optional
            # " Iterations: N." tail. The method label, the dispatch language
            # and the host version move with options.lang, and the iteration
            # count is a property of the run, not a model result.
            r'\w+ analysis \[[^\]]*\][^\n]*',
            # LINE's own diagnostics. warning() hides an identical repeat for
            # 60 seconds, so whether a banner or its suppression notice shows
            # up at all depends on how fast the run is rather than on the model.
            r'Warning \[\w+\]: [^\n]*',
            r'\[\w+\] Message cast more than once[^\n]*',
            # Notes on how the state space was built, cast by the native CTMC
            # path only.
            # The cutoff may be a per-(station,class) MATRIX, whose repr runs
            # over several lines, so this cannot stop at the first newline.
            r'(?s)CTMC solver using state space cutoff = .*?for open/mixed model\.',
            r'CTMC has \d+ SCCs[^\n]*',
            # State-space reduction note of the native CTMC path.
            r'Stochcomp: \d+ immediate states eliminated[^\n]*',
            # Command echoes (SolverJMT/SolverLDES/LQNS/...) quote randomly
            # named scratch directories.
            r'\w+ command: [^\n]*',
            # Python warning banners carry a source path and line number that
            # move with any edit to the emitting module.
            r'\S+:\d+: \w*Warning: [^\n]*',
            # Default object repr embeds a per-run memory address.
            r'at 0x[0-9a-fA-F]+',
            r'Execution time: \d+\.?\d*s?',
            r'Time elapsed: \d+\.?\d*s?',
            # tic/toc elapsed times the LQN examples print alongside results.
            r'\bT(?:no)?init\s*=\s*\d+\.?\d*',
            r'Runtime: \d+\.?\d*s?',
            r'Time: \d+\.?\d*s?',
            # A run of 10+ digits is a timestamp or a pid, but ONLY when it is
            # not part of a number: `\b\d{10,}\b` also matched the fractional
            # digits of a full-precision float, cutting 0.30769230769230765
            # down to 0. and reporting the cell's only result as 0.0.
            r'(?<![\d.])\d{10,}(?![\d.])',
        ]

        filtered_text = text
        for pattern in patterns_to_remove:
            filtered_text = re.sub(pattern, '', filtered_text, flags=re.IGNORECASE)

        return filtered_text

    def _compare_results(self, regression, current, filename):
        """Compare regression and current results, return list of differences."""
        differences = []

        for cell_key in regression:
            if cell_key not in current:
                differences.append(f"  {cell_key}: Missing in current results")
                continue

            regression_values = regression[cell_key]
            current_values = current[cell_key]

            if len(regression_values) != len(current_values):
                differences.append(
                    f"  {cell_key}: Different number of values "
                    f"(regression: {len(regression_values)}, current: {len(current_values)})"
                )
                continue

            for i, (regression_val, current_val) in enumerate(zip(regression_values, current_values)):
                abs_diff = abs(current_val - regression_val)
                rel_diff = abs_diff / abs(regression_val) if regression_val != 0 else abs_diff

                if abs_diff > TOLERANCE and rel_diff > RELATIVE_TOLERANCE:
                    differences.append(
                        f"  {cell_key}[{i}]: {current_val} vs {regression_val} "
                        f"(abs_diff: {abs_diff:.2e}, rel_diff: {rel_diff:.2e})"
                    )

        for cell_key in current:
            if cell_key not in regression:
                differences.append(f"  {cell_key}: Extra in current results")

        return differences

    def tearDown(self):
        """Print any baseline warnings at the end of each test."""
        if self.baseline_warnings:
            print("\n" + "="*60)
            print("BASELINE COMPARISON WARNINGS:")
            print("="*60)
            for warning in self.baseline_warnings:
                print(warning)
            print("="*60)

    @classmethod
    def examples_root(cls):
        working_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        return os.path.join(working_dir, cls.notebooks_dir)

    @classmethod
    def discover_notebooks(cls):
        """Every notebook under `roots`, as paths relative to notebooks_dir."""
        examples_dir = cls.examples_root()

        relative_paths = []
        for root in cls.roots:
            if root == '':
                # the examples directory itself, without descending
                pattern = os.path.join(examples_dir, '*.ipynb')
            else:
                root_path = os.path.join(examples_dir, root)
                if not os.path.isdir(root_path):
                    raise AssertionError(
                        f"{cls.__name__}: declared root '{root}' does not exist under "
                        f"{examples_dir}. Remove it or fix the name -- a root that "
                        f"silently vanishes takes its notebooks out of every suite.")
                pattern = os.path.join(root_path, '**', '*.ipynb')
            for filepath in glob.glob(pattern, recursive=True):
                if '.ipynb_checkpoints' in filepath:
                    continue
                relative_paths.append(os.path.relpath(filepath, examples_dir))

        relative_paths.sort()
        return relative_paths


def attach(cls):
    """Give `cls` one test method per notebook it discovers."""
    for notebook_path in cls.discover_notebooks():
        test_name = ('test_' + notebook_path.replace('/', '_')
                     .replace('.ipynb', '').replace('-', '_'))

        def _test_method(self, _path=notebook_path):
            self.run_notebook(_path)

        _test_method.__name__ = test_name
        _test_method.__doc__ = f"Test {notebook_path}"
        setattr(cls, test_name, _test_method)
    return cls
