"""Smoke tests for the root-level quickstart artifacts (mm1.py, mm1.ipynb).

These are the entry points advertised in python/README.md and CLAUDE.md. They live
outside examples/, so the notebook-driven suites (test_all_examples.py,
test_all_gettingstarted.py) never collect them. The checks below assert the
observable contract of the quickstart: a nonzero exit code is a failure, and the
avg table must reach stdout with M/M/1 values at rho = 0.5.
"""

import json
import os
import subprocess
import sys
import unittest

PYTHON_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# JMT is a simulator; the seed is fixed but the tolerance stays wide enough to
# absorb sampling noise on the default sample budget.
SIM_RTOL = 0.10

# M/M/1 with lambda = 1, mu = 2: rho = 0.5, QLen = 1, RespT = 1, Tput = 1.
EXPECTED = {'QLen': 1.0, 'Util': 0.5, 'RespT': 1.0, 'Tput': 1.0}


def parse_avg_table(stdout):
    """Parse the printed avg table into {station: {column: value}}."""
    rows = {}
    header = None
    for line in stdout.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'Station' and 'JobClass' in fields:
            header = fields[2:]
            continue
        if header is not None and len(fields) == len(header) + 2:
            try:
                values = [float(v) for v in fields[2:]]
            except ValueError:
                continue
            rows[fields[0]] = dict(zip(header, values))
    return rows


class TestQuickstartScript(unittest.TestCase):

    def run_script(self, name):
        proc = subprocess.run([sys.executable, name], cwd=PYTHON_ROOT,
                              capture_output=True, text=True, timeout=600)
        self.assertEqual(proc.returncode, 0,
                         f"{name} exited with {proc.returncode}\n{proc.stderr}")
        return proc.stdout

    def assert_mm1_table(self, rows, source_name, queue_name):
        self.assertIn(source_name, rows, f"no {source_name} row in the printed table")
        self.assertIn(queue_name, rows, f"no {queue_name} row in the printed table")
        for column, expected in EXPECTED.items():
            self.assertIn(column, rows[queue_name], f"column {column} missing")
            self.assertAlmostEqual(rows[queue_name][column], expected,
                                   delta=SIM_RTOL * expected,
                                   msg=f"{queue_name}.{column} off the M/M/1 value")

    def test_mm1_script_prints_avg_table(self):
        stdout = self.run_script('mm1.py')
        rows = parse_avg_table(stdout)
        self.assertTrue(rows, f"mm1.py printed no avg table:\n{stdout}")
        self.assert_mm1_table(rows, 'Source', 'Queue')

    def test_mm1_notebook_prints_avg_table(self):
        try:
            import nbformat
            from nbconvert.preprocessors import ExecutePreprocessor
        except ImportError:
            self.skipTest('nbformat/nbconvert unavailable')

        path = os.path.join(PYTHON_ROOT, 'mm1.ipynb')
        with open(path) as f:
            nb = nbformat.read(f, as_version=4)
        # THE KERNEL BOOT IS A SECOND DEADLINE, and nbclient defaults it to 60 s.
        # It runs before the first cell, so `timeout` never covers it, and under
        # `-n auto` on a loaded host 60 s is a contention detector rather than a
        # hang detector -- see notebook_suite.KERNEL_STARTUP, which carries the
        # measurement and lost two notebooks to this on 2026-09-09. Spelled out
        # here rather than imported: importing that module installs a kernelspec
        # under an flock, and this test drives the plain `python3` kernel.
        ExecutePreprocessor(timeout=600, kernel_name='python3',
                            startup_timeout=480).preprocess(
            nb, {'metadata': {'path': PYTHON_ROOT}})

        stdout = ''
        for cell in nb.cells:
            for output in cell.get('outputs', []):
                self.assertNotEqual(output.output_type, 'error',
                                    f"mm1.ipynb raised {output.get('ename')}: {output.get('evalue')}")
                if output.output_type == 'stream':
                    stdout += output.get('text', '')
                elif output.output_type == 'execute_result':
                    stdout += output.get('data', {}).get('text/plain', '')

        rows = parse_avg_table(stdout)
        self.assertTrue(rows, f"mm1.ipynb printed no avg table:\n{stdout}")
        # Same node names as the script: mm1.ipynb is generated from mm1.py by
        # tools/gen_notebooks.py, so the two cannot name their nodes differently.
        self.assert_mm1_table(rows, 'Source', 'Queue')

    def test_readme_documents_the_printed_columns(self):
        """The README example block must list the columns the script really prints."""
        with open(os.path.join(PYTHON_ROOT, 'README.md')) as f:
            readme = f.read()
        stdout = self.run_script('mm1.py')
        rows = parse_avg_table(stdout)
        self.assertTrue(rows, f"mm1.py printed no avg table:\n{stdout}")
        documented = parse_avg_table(readme)
        self.assertTrue(documented, 'README.md has no avg table example block')
        printed_columns = set(next(iter(rows.values())))
        for station, columns in documented.items():
            self.assertIn(station, rows, f"README documents an absent station {station}")
            self.assertEqual(set(columns), printed_columns,
                             'README example block columns are stale')


if __name__ == '__main__':
    unittest.main()
