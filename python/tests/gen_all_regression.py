
import json
import os
import shutil
import sys
import traceback
from pathlib import Path
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import re
import numpy as np
import warnings

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# The single source of truth for which notebooks exist and where their baseline
# lives. This file used to carry its own hardcoded subfolder list, which drifted
# out of step with the suites' own list -- a notebook could be executed by a
# suite and yet have no baseline anyone could generate for it.
from notebook_suite import NotebookSuite

SOURCE_DIR = 'examples'


class RegressionGenerator:
    def __init__(self):
        self.regression_dir = 'tests/regression'
        self.working_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

        regression_path = os.path.join(self.working_dir, self.regression_dir)
        os.makedirs(regression_path, exist_ok=True)

    def extract_cell_content(self, nb):
        """Extract comprehensive cell content from notebook including source, outputs, and metadata."""
        results = {}

        for cell_idx, cell in enumerate(nb.cells):
            cell_data = {
                'cell_type': cell.cell_type,
                'source': cell.get('source', ''),
                'metadata': cell.get('metadata', {}),
                'outputs': [],
                'execution_count': cell.get('execution_count', None)
            }

            if cell.cell_type == 'code':
                for output_idx, output in enumerate(cell.get('outputs', [])):
                    output_data = {
                        'output_type': output.output_type,
                        'execution_count': output.get('execution_count', None),
                        'metadata': output.get('metadata', {}),
                        'content': {}
                    }

                    if output.output_type in ['execute_result', 'display_data']:
                        if 'data' in output:
                            for data_type, data_content in output['data'].items():
                                if data_type == 'text/plain':
                                    output_data['content'][data_type] = data_content
                                    numerical_values = self._extract_numbers_from_text(data_content)
                                    if numerical_values:
                                        output_data['content']['numerical_values'] = numerical_values
                                elif data_type in ['text/html', 'application/json', 'image/png', 'image/jpeg']:
                                    if data_type.startswith('image/'):
                                        output_data['content'][data_type] = f"<image_data_length:{len(str(data_content))}>"
                                    else:
                                        output_data['content'][data_type] = data_content

                    elif output.output_type == 'stream':
                        stream_name = output.get('name', 'unknown')
                        text_content = output.get('text', '')
                        output_data['content']['stream_name'] = stream_name
                        output_data['content']['text'] = text_content

                        numerical_values = self._extract_numbers_from_text(text_content)
                        if numerical_values:
                            output_data['content']['numerical_values'] = numerical_values

                    elif output.output_type == 'error':
                        output_data['content']['ename'] = output.get('ename', '')
                        output_data['content']['evalue'] = output.get('evalue', '')
                        output_data['content']['traceback'] = output.get('traceback', [])

                    cell_data['outputs'].append(output_data)

            if cell_data['source'] or cell_data['outputs']:
                results[f'cell_{cell_idx}'] = cell_data

        return results

    def _extract_numbers_from_text(self, text):
        """Extract numerical values from text output using regex."""
        if not isinstance(text, str):
            if isinstance(text, (list, tuple)):
                text = ' '.join(str(item) for item in text)
            else:
                return []

        number_pattern = r'[-+]?(?:(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)'

        matches = re.findall(number_pattern, text)

        numbers = []
        for match in matches:
            try:
                num = float(match)
                if np.isfinite(num):
                    numbers.append(num)
            except (ValueError, OverflowError):
                continue

        return numbers

    def extract_numerical_summary(self, cell_content):
        """Extract a summary of numerical values for backwards compatibility."""
        numerical_summary = {}

        for cell_idx, cell_data in cell_content.items():
            cell_numbers = []

            for output in cell_data.get('outputs', []):
                if 'numerical_values' in output.get('content', {}):
                    cell_numbers.extend(output['content']['numerical_values'])

            source = cell_data.get('source', '')
            if source:
                source_numbers = self._extract_numbers_from_text(source)
                if len(source_numbers) <= 10:
                    cell_numbers.extend(source_numbers)

            if cell_numbers:
                numerical_summary[cell_idx] = cell_numbers

        return numerical_summary

    def generate_regression(self, notebook_path, source_dir):
        """Generate regression results for a single notebook."""
        print(f'Generating regression data for: {source_dir}/{notebook_path}')

        os.chdir(self.working_dir)

        full_notebook_path = os.path.join(source_dir, notebook_path)

        try:
            with open(full_notebook_path) as f:
                nb = nbformat.read(f, as_version=4)

            ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
            ep.preprocess(nb, {'metadata': {'path': os.path.dirname(full_notebook_path)}})

            full_results = self.extract_cell_content(nb)

            results = {
                'cell_content': full_results,
                'numerical_summary': self.extract_numerical_summary(full_results),
                'metadata': {
                    'notebook_path': notebook_path,
                    'source_dir': source_dir,
                    'total_cells': len(nb.cells),
                    'code_cells': len([cell for cell in nb.cells if cell.cell_type == 'code']),
                    'generation_info': {
                        'script_version': '2.0',
                        'extraction_method': 'comprehensive'
                    }
                }
            }

            notebook_dir = os.path.dirname(notebook_path)
            notebook_name = os.path.basename(notebook_path).replace('.ipynb', '_regression.json')

            regression_subfolder = os.path.join(self.working_dir, self.regression_dir, source_dir, notebook_dir)
            os.makedirs(regression_subfolder, exist_ok=True)

            regression_path = os.path.join(regression_subfolder, notebook_name)

            with open(regression_path, 'w') as f:
                json.dump(results, f, indent=2, default=self._json_serializer)

            print(f'✓ Regression data saved: {regression_path}')
            return True

        except Exception as e:
            print(f'✗ Failed to generate regression data for {source_dir}/{notebook_path}: {e}')
            traceback.print_exc()
            return False

        finally:
            if os.path.exists("_trial_temp"):
                shutil.rmtree("_trial_temp")

    def _json_serializer(self, obj):
        """Custom JSON serializer for numpy types and other objects."""
        if isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif hasattr(obj, 'tolist'):
            return obj.tolist()
        elif hasattr(obj, '__dict__'):
            return str(obj)
        raise TypeError(f"Object {obj} of type {type(obj)} is not JSON serializable")

    def generate_all_regression(self):
        """Generate regression results for every notebook the suites execute."""
        print("Generating regression data for all Python notebooks...")
        print(f"Working directory: {self.working_dir}")
        print(f"Source directory: {SOURCE_DIR}")
        print(f"Regression data will be saved to: {os.path.join(self.working_dir, self.regression_dir)}")

        notebook_files = all_notebooks()
        print(f"Found {len(notebook_files)} notebook files in {SOURCE_DIR}/")

        total_successful = 0
        total_failed = 0

        for notebook_file in notebook_files:
            if self.generate_regression(notebook_file, SOURCE_DIR):
                total_successful += 1
            else:
                total_failed += 1

        print(f"\n{'='*60}")
        print(f"OVERALL REGRESSION GENERATION SUMMARY")
        print(f"{'='*60}")
        print(f"  Total successful: {total_successful}")
        print(f"  Total failed: {total_failed}")
        print(f"  Grand total: {total_successful + total_failed}")

        if total_failed > 0:
            print(f"\nWarning: {total_failed} notebooks failed to generate regression data.")
            return False

        return True

def all_notebooks():
    """Every notebook the per-area suites execute, relative to `examples/`."""

    class _All(NotebookSuite):
        roots = ['basic', 'advanced', 'gettingstarted', 'gallery',
                 'discrete', 'inference', 'opt', 'solvers', '']

    return _All.discover_notebooks()


def _split_source_dir(path):
    """Split a notebook path into the (source_dir, relative path) pair used above.

    Accepts either a path rooted at python/ ('examples/basic/x/y.ipynb') or one
    already relative to `examples/` ('basic/x/y.ipynb'). Anything the suites do
    not execute has no baseline to write and is rejected rather than written to
    a path nothing reads.
    """
    path = path.replace(os.sep, '/').lstrip('./')
    if path.startswith(SOURCE_DIR + '/'):
        path = path[len(SOURCE_DIR) + 1:]
    if path in all_notebooks():
        return SOURCE_DIR, path
    return None, path


def main():
    """Main function to run regression generation.

    With notebook paths on the command line only those baselines are rewritten.
    Rebasing every notebook at once is almost never right: a baseline that moved
    because a solver was fixed and one that moved because a solver broke look
    identical here, so the whole set has to be adjudicated notebook by notebook
    and only the adjudicated ones regenerated.
    """
    print("=" * 60)
    print("LINE Solver Python Notebooks - Regression Generator")
    print("=" * 60)

    generator = RegressionGenerator()

    targets = sys.argv[1:]
    if targets:
        success = True
        for target in targets:
            source_dir, rel_path = _split_source_dir(target)
            if source_dir is None:
                print(f'✗ {target} is not a notebook any example suite executes; '
                      f'expected a path under examples/')
                success = False
                continue
            success &= generator.generate_regression(rel_path, source_dir)
        sys.exit(0 if success else 1)

    success = generator.generate_all_regression()

    if success:
        print("\n✓ All regression data generated successfully!")
        print("You can now run test_all_examples.py to verify notebooks against this regression data.")
        print("\nRegression data structure:")
        print("  tests/regression/examples/<area>/<name>_regression.json")
    else:
        print("\n✗ Some regression data failed to generate.")
        print("Check the error messages above and fix any issues before running tests.")
        sys.exit(1)

if __name__ == '__main__':
    main()