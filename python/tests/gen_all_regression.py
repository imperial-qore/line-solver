
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

class RegressionGenerator:
    def __init__(self):
        self.notebook_dirs = ['examples/gettingstarted', 'examples/basic', 'examples/advanced']
        self.regression_dir = 'tests/regression'
        self.working_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

        self.enabled_examples_subfolders = [
            'cacheModel',
            'cdfRespT',
            'classSwitching',
            'closedQN',
            'cyclicPolling',
            'forkJoin',
            'initState',
            'layeredCQ',
            'layeredModel',
            'loadDependent',
            'misc',
            'mixedQN',
            'openQN',
            'prioModel',
            'randomEnv',
            'rewardModel',
            'stateDepRouting',
            'stateProbabilities',
            'stochPetriNet',
            'switchoverTimes',
        ]

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
        """Generate regression results for all notebooks in advanced and gettingstarted."""
        print("Generating regression data for all Python notebooks...")
        print(f"Working directory: {self.working_dir}")
        print(f"Source directories: {', '.join(self.notebook_dirs)}")
        print(f"Regression data will be saved to: {os.path.join(self.working_dir, self.regression_dir)}")

        total_successful = 0
        total_failed = 0

        for source_dir in self.notebook_dirs:
            print(f"\n{'='*50}")
            print(f"Processing {source_dir}/ directory")
            print(f"{'='*50}")

            source_path = os.path.join(self.working_dir, source_dir)
            if not os.path.exists(source_path):
                print(f"Warning: Directory {source_path} does not exist, skipping...")
                continue

            notebook_files = []

            if source_dir in ['examples/advanced', 'examples/basic']:
                for subfolder in self.enabled_examples_subfolders:
                    subfolder_path = os.path.join(source_path, subfolder)
                    if not os.path.exists(subfolder_path):
                        continue

                    for file in os.listdir(subfolder_path):
                        if file.endswith('.ipynb'):
                            rel_path = os.path.join(subfolder, file)
                            notebook_files.append(rel_path)
            else:
                for root, dirs, files in os.walk(source_path):
                    for file in files:
                        if file.endswith('.ipynb'):
                            rel_path = os.path.relpath(os.path.join(root, file), source_path)
                            notebook_files.append(rel_path)

            notebook_files.sort()

            print(f"Found {len(notebook_files)} notebook files in {source_dir}/")

            successful = 0
            failed = 0

            for notebook_file in notebook_files:
                if self.generate_regression(notebook_file, source_dir):
                    successful += 1
                else:
                    failed += 1

            print(f"\n{source_dir}/ results:")
            print(f"  Successful: {successful}")
            print(f"  Failed: {failed}")
            print(f"  Total: {len(notebook_files)}")

            total_successful += successful
            total_failed += failed

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

def main():
    """Main function to run regression generation."""
    print("=" * 60)
    print("LINE Solver Python Notebooks - Regression Generator")
    print("=" * 60)

    generator = RegressionGenerator()
    success = generator.generate_all_regression()

    if success:
        print("\n✓ All regression data generated successfully!")
        print("You can now run test_all_examples.py to verify notebooks against this regression data.")
        print("\nRegression data structure:")
        print("  tests/regression/advanced/     - Examples regression data")
        print("  tests/regression/gettingstarted/ - Getting started regression data")
    else:
        print("\n✗ Some regression data failed to generate.")
        print("Check the error messages above and fix any issues before running tests.")
        sys.exit(1)

if __name__ == '__main__':
    main()