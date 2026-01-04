import shutil
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import os
import sys
import unittest
import json
import warnings
import re
import numpy as np
from pathlib import Path
import glob
import traceback
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Install a kernel for the current Python interpreter to ensure notebooks run with the same environment
KERNEL_NAME = 'line_solver_test_kernel'
try:
    subprocess.run(
        [sys.executable, '-m', 'ipykernel', 'install', '--user', '--name', KERNEL_NAME, '--display-name', 'LINE Test Kernel'],
        check=True, capture_output=True
    )
except subprocess.CalledProcessError:
    # Fall back to default python3 kernel if installation fails
    KERNEL_NAME = 'python3'

TOLERANCE = 1e-6
RELATIVE_TOLERANCE = 1e-5

class TestExamples(unittest.TestCase):

    notebooks_dir = 'examples'
    baselines_dir = 'tests/regression'

    enabled_subfolders = [
        'advanced/cdfRespT',
        'advanced/cyclicPolling',
        'advanced/initState',
        'advanced/layeredCQ',
        'advanced/loadDependent',
        'advanced/randomEnv',
        'advanced/stateDepRouting',
        'advanced/stateProbabilities',
        'advanced/switchoverTimes',
    ]
#    'advanced/misc',
#    'advanced/rewardModel',

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
            ep = ExecutePreprocessor(timeout=600, kernel_name=KERNEL_NAME)
            execution_error = None
            try:
                ep.preprocess(nb, {'metadata': {'path': os.path.dirname(filepath)}})
                print(f'✓ Successfully executed notebook: {filename}', flush=True)
            except Exception as e:
                execution_error = e
                print(f"\n{'='*60}", flush=True)
                print(f"NOTEBOOK EXECUTION FAILED: {filename}", flush=True)
                print(f"{'='*60}", flush=True)
                print(f"Error Type: {type(e).__name__}", flush=True)
                print(f"Error Message: {str(e)}", flush=True)
                print(f"\nFull Traceback:", flush=True)
                print(traceback.format_exc(), flush=True)
                print(f"{'='*60}", flush=True)

                raise

            finally:
                pass

        for cell_idx, cell in enumerate(nb.cells):
            if cell.cell_type == 'code':
                for output in cell.get('outputs', []):
                    if output.output_type == 'error':
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
        regression_path = os.path.join(self.working_dir, self.baselines_dir, self.notebooks_dir, notebook_dir, notebook_name)

        if not os.path.exists(regression_path):
            warning_msg = f"No regression data found for {filename}. Run gen_all_regression.py to generate regression data."
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

        current_results = self._extract_numerical_outputs(nb)

        differences = self._compare_results(regression_results, current_results, filename)

        if differences:
            warning_msg = f"Numerical differences detected in {filename}:\n" + "\n".join(differences)
            warnings.warn(warning_msg, UserWarning)
            self.baseline_warnings.append(warning_msg)


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
        """Extract numerical values from text output using regex."""
        if not isinstance(text, str):
            return []

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
        """Print any baseline warnings and prompt generation summary at the end of each test."""
        if self.baseline_warnings:
            print("\n" + "="*60)
            print("BASELINE COMPARISON WARNINGS:")
            print("="*60)
            for warning in self.baseline_warnings:
                print(warning)
            print("="*60)


    @classmethod
    def discover_notebooks(cls):
        """Discover notebooks in enabled subfolders only."""
        working_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        examples_dir = os.path.join(working_dir, cls.notebooks_dir)

        relative_paths = []

        for subfolder in cls.enabled_subfolders:
            subfolder_path = os.path.join(examples_dir, subfolder)
            if not os.path.exists(subfolder_path):
                print(f"Warning: Enabled subfolder {subfolder} does not exist, skipping...")
                continue

            pattern = os.path.join(subfolder_path, '*.ipynb')
            notebook_files = glob.glob(pattern)

            for filepath in notebook_files:
                relative_path = os.path.relpath(filepath, examples_dir)
                relative_paths.append(relative_path)

        relative_paths.sort()
        return relative_paths

def create_test_method(notebook_path):
    """Create a test method for a specific notebook."""
    def _test_method(self):
        self.run_notebook(notebook_path)
    return _test_method

for notebook_path in TestExamples.discover_notebooks():
    test_name = 'test_' + notebook_path.replace('/', '_').replace('.ipynb', '').replace('-', '_')

    test_method = create_test_method(notebook_path)
    test_method.__name__ = test_name
    test_method.__doc__ = f"Test {notebook_path}"
    setattr(TestExamples, test_name, test_method)

# Clean up module-level variables to prevent pytest from picking them up as tests
del notebook_path, test_name, test_method

if __name__ == '__main__':
    print(f"Running tests with tolerance settings:")
    print(f"  Absolute tolerance: {TOLERANCE}")
    print(f"  Relative tolerance: {RELATIVE_TOLERANCE}")
    print(f"  Regression directory: tests/regression/")
    print("="*60)
    unittest.main(verbosity=2)
