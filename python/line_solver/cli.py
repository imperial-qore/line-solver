"""
LINE Command-Line Interface

Provides command-line and WebSocket server interface for LINE solver functionality.
All solver requests are routed through the LINE class (SolverAuto) which automatically
selects the best solver or uses a specified solver/method combination.

Supports loading models from multiple formats (JSIM, LQN, MAT, pickle),
running analysis with various solvers, and outputting results in multiple formats.

Usage:
    line --help                              Show help
    line --file model.lqn --solver mva       Solve LQN model with MVA
    line --file model.jsimg --solver mva.lin Use MVA with linearizer method
    line --port 5863                         Start WebSocket server
    cat model.jsimg | line --input jsimg     Solve JSIM from stdin
"""

import sys
import os
import argparse
import json
import csv
import pickle
import tempfile
import asyncio
from pathlib import Path
from io import StringIO

# Import LINE solver components
from line_solver import LINE, VerboseLevel, GlobalConstants
from line_solver.io import M2M
from line_solver.layered import LayeredNetwork


# Valid solver prefixes (for help text)
VALID_SOLVERS = ['auto', 'line', 'sim', 'exact', 'fast', 'accurate',
                 'mva', 'ctmc', 'fluid', 'jmt', 'nc', 'ssa', 'ln', 'lqns',
                 'mam', 'des', 'env', 'qns', 'posterior']

# Example solver.method combinations
SOLVER_METHOD_EXAMPLES = [
    'auto', 'sim', 'exact', 'fast', 'accurate',
    'mva', 'mva.lin', 'mva.exact', 'mva.amva',
    'nc', 'nc.exact', 'nc.ls',
    'ctmc', 'ctmc.gpu',
    'fluid', 'fluid.statedep',
    'jmt', 'jmt.jsim', 'jmt.jmva',
]

# Solver compatibility with input formats
JSIM_COMPATIBLE_SOLVERS = ['ctmc', 'fluid', 'jmt', 'mva', 'nc', 'ssa', 'mam', 'des', 'auto']
LQN_COMPATIBLE_SOLVERS = ['ln', 'lqns', 'mva', 'nc', 'env']
JSIM_FORMATS = ['jsim', 'jsimg', 'jsimw']
LQN_FORMATS = ['lqnx', 'xml']
MAT_FORMATS = ['mat']
SUPPORTED_FORMATS = JSIM_FORMATS + LQN_FORMATS + MAT_FORMATS + ['pkl']

# Output formats
OUTPUT_FORMATS = ['readable', 'json', 'csv', 'pickle', 'mat']

# All valid analysis types
VALID_ANALYSIS_TYPES = [
    # Basic
    'all', 'avg', 'sys', 'stage', 'chain', 'node', 'nodechain',
    # Distribution
    'cdf-respt', 'cdf-passt', 'perct-respt',
    # Transient
    'tran-avg', 'tran-cdf-respt', 'tran-cdf-passt',
    # Probability
    'prob', 'prob-aggr', 'prob-marg', 'prob-sys', 'prob-sys-aggr',
    # Sampling
    'sample', 'sample-aggr', 'sample-sys', 'sample-sys-aggr',
    # Reward
    'reward', 'reward-steady', 'reward-value'
]

# Analysis types that require specific solvers
ANALYSIS_SOLVER_COMPAT = {
    'sample': ['ssa'],
    'sample-aggr': ['ssa'],
    'sample-sys': ['ssa'],
    'sample-sys-aggr': ['ssa'],
    'reward': ['ctmc'],
    'reward-steady': ['ctmc'],
    'reward-value': ['ctmc'],
    'perct-respt': ['mam'],
    'prob': ['ctmc', 'ssa'],
    'prob-aggr': ['ctmc', 'ssa'],
    'prob-marg': ['ctmc', 'ssa'],
    'prob-sys': ['ctmc', 'ssa'],
    'prob-sys-aggr': ['ctmc', 'ssa'],
}

# Analysis types that require node index
ANALYSIS_REQUIRES_NODE = ['prob', 'prob-aggr', 'prob-marg', 'sample', 'sample-aggr']

# Analysis types that require class index
ANALYSIS_REQUIRES_CLASS = ['prob-marg']

# Valid built-in reward names
VALID_REWARD_NAMES = ['QLen', 'Tput', 'Util', 'RespT', 'WaitT', 'ArvR', 'ResidT']


def create_parser():
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog='line',
        description='LINE - Command-Line Interface for Queueing Network Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with file input
  line -f model.lqn -s mva -a avg

  # Using solver.method syntax
  line -f model.jsimg -s mva.lin
  line -f model.jsimg -s nc.exact

  # High-level method selection
  line -f model.jsimg -s sim       # Best simulator
  line -f model.jsimg -s exact     # Best exact method
  line -f model.jsimg -s fast      # Fast approximate
  line -f model.jsimg -s accurate  # Accurate approximate

  # Server mode
  line -p 8080

  # LQN model with layered network solver
  line -f model.lqnx -s ln

Solver.Method Syntax:
  Solvers: mva, nc, ctmc, fluid, jmt, ssa, mam, des, ln, lqns, env
  Methods vary by solver, e.g.:
    mva.lin, mva.exact, mva.amva
    nc.exact, nc.ls, nc.comom
    ctmc.gpu
    fluid.statedep, fluid.closing
    jmt.jsim, jmt.jmva
        """
    )

    # Port/server mode
    parser.add_argument(
        '-p', '--port',
        type=int,
        default=None,
        metavar='PORT',
        help='Run LINE solver in server mode on specified port (default: 5863)'
    )

    # Input file
    parser.add_argument(
        '-f', '--file',
        type=str,
        default=None,
        metavar='PATH',
        help='Specify input model file (default: read from stdin)'
    )

    # Input format
    parser.add_argument(
        '-i', '--input',
        type=str,
        choices=SUPPORTED_FORMATS,
        default=None,
        metavar='FORMAT',
        help=f'Input file format: {", ".join(SUPPORTED_FORMATS)} (default: auto-detect from extension)'
    )

    # Output format
    parser.add_argument(
        '-o', '--output',
        type=str,
        choices=OUTPUT_FORMATS,
        default='readable',
        metavar='FORMAT',
        help=f'Output format: {", ".join(OUTPUT_FORMATS)} (default: readable)'
    )

    # Solver selection
    parser.add_argument(
        '-s', '--solver',
        type=str,
        default='auto',
        metavar='SOLVER',
        help='''Solver or solver.method (default: auto). Options:
High-level: auto, sim, exact, fast, accurate
Solvers: mva, nc, ctmc, fluid, jmt, ssa, mam, des, ln, lqns
With method: mva.lin, nc.exact, ctmc.gpu, etc.'''
    )

    # Analysis type (supports comma-separated values)
    parser.add_argument(
        '-a', '--analysis',
        type=str,
        default='all',
        metavar='TYPE',
        help='''Analysis type(s), comma-separated for multiple. Options:
Basic: all, avg, sys, stage, chain, node, nodechain
CDF: cdf-respt, cdf-passt, perct-respt
Transient: tran-avg, tran-cdf-respt, tran-cdf-passt
Probability: prob, prob-aggr, prob-marg, prob-sys, prob-sys-aggr
Sampling: sample, sample-aggr, sample-sys, sample-sys-aggr
Reward: reward, reward-steady, reward-value (default: all)'''
    )

    # Node index for prob/sample analysis
    parser.add_argument(
        '-n', '--node',
        type=int,
        default=None,
        metavar='INDEX',
        help='Node index for prob/sample analysis (0-based)'
    )

    # Class index for prob-marg analysis
    parser.add_argument(
        '-c', '--class-idx',
        type=int,
        default=None,
        metavar='INDEX',
        help='Job class index for prob-marg analysis (0-based)'
    )

    # State vector for prob analysis
    parser.add_argument(
        '--state',
        type=str,
        default=None,
        metavar='VALUES',
        help='State vector for prob analysis (comma-separated integers)'
    )

    # Number of events for sample analysis
    parser.add_argument(
        '--events',
        type=int,
        default=1000,
        metavar='NUMBER',
        help='Number of events for sample analysis (default: 1000)'
    )

    # Percentiles for perct-respt
    parser.add_argument(
        '--percentiles',
        type=str,
        default='50,90,95,99',
        metavar='VALUES',
        help='Percentile values for perct-respt (comma-separated, default: 50,90,95,99)'
    )

    # Reward name for reward-value
    parser.add_argument(
        '--reward-name',
        type=str,
        default=None,
        choices=VALID_REWARD_NAMES,
        metavar='NAME',
        help=f'Built-in reward name for reward-value analysis: {", ".join(VALID_REWARD_NAMES)}'
    )

    # Random seed
    parser.add_argument(
        '-d', '--seed',
        type=int,
        default=None,
        metavar='SEED',
        help='Random seed for stochastic solvers (JMT, SSA, DES)'
    )

    # Verbosity
    parser.add_argument(
        '-v', '--verbosity',
        type=str,
        choices=['normal', 'silent'],
        default='normal',
        metavar='LEVEL',
        help='Verbosity level: normal (default) or silent'
    )

    # Docker mode
    parser.add_argument(
        '--docker-host',
        type=str,
        default='127.0.0.1',
        metavar='HOST',
        help='Docker server host (default: 127.0.0.1)'
    )

    parser.add_argument(
        '--docker-port',
        type=str,
        default='5863',
        metavar='PORT',
        help='Docker server port (default: 5863)'
    )

    # Version
    try:
        import line_solver
        version = getattr(line_solver, '__version__', '3.0.3')
    except:
        version = '3.0.3'

    parser.add_argument(
        '-V', '--version',
        action='version',
        version=f'LINE Solver {version}'
    )

    return parser


def detect_input_format(filename):
    """Auto-detect input format from filename extension."""
    if filename is None:
        return None

    _, ext = os.path.splitext(filename.lower())
    ext = ext.lstrip('.')

    if ext in SUPPORTED_FORMATS:
        return ext

    raise ValueError(f"Cannot auto-detect format from extension: {ext}")


def load_model(filename, input_format, verbose=False):
    """Load model from file, detecting format if not specified."""
    if filename is None:
        raise ValueError("Filename required for model loading")

    # Auto-detect format if not specified
    if input_format is None:
        input_format = detect_input_format(filename)
        if verbose:
            print(f"Auto-detected format: {input_format}")

    # Check file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Model file not found: {filename}")

    # Load based on format
    if input_format in JSIM_FORMATS:
        m2m = M2M()
        model = m2m.JSIM2LINE(filename)
        if verbose:
            print(f"Loaded JSIM model from: {filename}")
        return model

    elif input_format in LQN_FORMATS:
        model = LayeredNetwork.load(filename, verbose)
        if verbose:
            print(f"Loaded LQN model from: {filename}")
        return model

    elif input_format == 'mat':
        m2m = M2M()
        model = m2m.MAT2LINE(filename)
        if verbose:
            print(f"Loaded MATLAB model from: {filename}")
        return model

    elif input_format == 'pkl':
        with open(filename, 'rb') as f:
            model = pickle.load(f)
        if verbose:
            print(f"Loaded pickled model from: {filename}")
        return model

    else:
        raise ValueError(f"Unsupported format: {input_format}")


def load_from_stdin(input_format=None, verbose=False):
    """Load model from stdin and save to temp file."""
    content = sys.stdin.read()

    if not content.strip():
        raise ValueError("No model content provided on stdin")

    # Determine format if not specified
    if input_format is None:
        input_format = 'jsim'  # Default to JSIM
        if verbose:
            print(f"No format specified, defaulting to: {input_format}")

    # Create temporary file with appropriate extension
    suffix = f'.{input_format}'
    with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as f:
        f.write(content)
        temp_filename = f.name

    try:
        model = load_model(temp_filename, input_format, verbose)
        return model
    finally:
        # Clean up temp file
        try:
            os.unlink(temp_filename)
        except:
            pass


def parse_solver_method(solver_str):
    """Parse 'solver.method' syntax into (solver, method) tuple.

    Examples:
        'auto' -> ('auto', 'default')
        'mva' -> ('mva', 'default')
        'mva.lin' -> ('mva', 'lin')
        'nc.exact' -> ('nc', 'exact')
    """
    if '.' in solver_str:
        parts = solver_str.split('.', 1)
        return parts[0], parts[1]
    return solver_str, 'default'


def validate_solver_compatibility(input_format, solver_str):
    """Validate that solver is compatible with input format."""
    # Extract base solver from solver.method syntax
    solver, _ = parse_solver_method(solver_str)

    # High-level methods are always valid (LINE handles routing)
    if solver in ['auto', 'line', 'sim', 'exact', 'fast', 'accurate']:
        return

    if input_format in JSIM_FORMATS:
        if solver not in JSIM_COMPATIBLE_SOLVERS:
            raise ValueError(
                f"Solver '{solver}' is not compatible with JSIM format. "
                f"Valid solvers: {', '.join(JSIM_COMPATIBLE_SOLVERS)}"
            )
    elif input_format in LQN_FORMATS:
        if solver not in LQN_COMPATIBLE_SOLVERS:
            raise ValueError(
                f"Solver '{solver}' is not compatible with LQN format. "
                f"Valid solvers: {', '.join(LQN_COMPATIBLE_SOLVERS)}"
            )


def solve_model(model, solver_str, seed=None, verbose=False):
    """Instantiate solver via LINE (SolverAuto) and run analysis.

    All solver requests are routed through LINE which handles:
    - Automatic solver selection (auto, sim, exact, fast, accurate)
    - Specific solver targeting (mva, nc, ctmc, jmt, etc.)
    - Solver.method syntax (mva.lin, nc.exact, ctmc.gpu, etc.)

    Args:
        model: The queueing network model to solve
        solver_str: Solver specification (e.g., 'auto', 'mva', 'mva.lin')
        seed: Random seed for stochastic solvers
        verbose: Enable verbose output

    Returns:
        Configured solver instance ready for analysis
    """
    if verbose:
        print(f"Creating solver with method='{solver_str}'...")

    # Build options for LINE solver
    options = {'verbose': verbose}

    if seed is not None:
        options['seed'] = seed

    # Route all requests through LINE (SolverAuto)
    # LINE handles solver selection based on the method parameter
    solver = LINE(model, method=solver_str, **options)

    if verbose:
        print("Running analysis...")

    return solver


def validate_analysis_types(analysis_str):
    """Validate comma-separated analysis types."""
    types = [t.strip() for t in analysis_str.split(',')]
    for t in types:
        if t not in VALID_ANALYSIS_TYPES:
            raise ValueError(f"Invalid analysis type: '{t}'. Valid types: {', '.join(VALID_ANALYSIS_TYPES)}")
    return types


def validate_analysis_solver_compat(analysis_types, solver_str):
    """Validate that analysis types are compatible with solver."""
    # Extract base solver from solver.method syntax
    solver_name, _ = parse_solver_method(solver_str)

    # Skip validation for high-level methods (LINE will handle routing)
    if solver_name in ['auto', 'line', 'sim', 'exact', 'fast', 'accurate']:
        return

    for analysis_type in analysis_types:
        required_solvers = ANALYSIS_SOLVER_COMPAT.get(analysis_type)
        if required_solvers and solver_name not in required_solvers:
            raise ValueError(
                f"Analysis type '{analysis_type}' requires solver: {' or '.join(required_solvers)}, "
                f"but '{solver_name}' was specified."
            )


def validate_analysis_params(analysis_types, node_idx, class_idx, reward_name):
    """Validate that required parameters are provided for analysis types."""
    for analysis_type in analysis_types:
        if analysis_type in ANALYSIS_REQUIRES_NODE and node_idx is None:
            raise ValueError(f"Analysis type '{analysis_type}' requires --node parameter.")
        if analysis_type in ANALYSIS_REQUIRES_CLASS and class_idx is None:
            raise ValueError(f"Analysis type '{analysis_type}' requires --class-idx parameter.")
        if analysis_type == 'reward-value' and not reward_name:
            raise ValueError("Analysis type 'reward-value' requires --reward-name parameter.")


def parse_state(state_str):
    """Parse state vector from comma-separated string."""
    if state_str is None:
        return None
    return [int(x.strip()) for x in state_str.split(',')]


def parse_percentiles(percentiles_str):
    """Parse percentiles from comma-separated string."""
    return [float(x.strip()) for x in percentiles_str.split(',')]


def execute_single_analysis(solver, model, analysis_type, node_idx, class_idx, state, num_events, percentiles, reward_name):
    """Execute a single analysis type and return results."""
    try:
        # Basic analysis types
        if analysis_type == 'all':
            results = {}
            try:
                results['avg'] = solver.avg_table()
            except:
                pass
            try:
                if hasattr(solver, 'avg_sys_table'):
                    results['sys'] = solver.avg_sys_table()
            except:
                pass
            return results

        elif analysis_type == 'avg':
            return solver.avg_table()

        elif analysis_type == 'sys':
            if hasattr(solver, 'avg_sys_table'):
                return solver.avg_sys_table()
            return None

        elif analysis_type == 'stage':
            if hasattr(solver, 'get_stage_table'):
                return solver.get_stage_table()
            return None

        elif analysis_type == 'chain':
            if hasattr(solver, 'get_avg_chain_table'):
                return solver.get_avg_chain_table()
            return None

        elif analysis_type == 'node':
            if hasattr(solver, 'get_avg_node_table'):
                return solver.get_avg_node_table()
            return None

        elif analysis_type == 'nodechain':
            if hasattr(solver, 'get_avg_node_chain_table'):
                return solver.get_avg_node_chain_table()
            return None

        # Distribution analysis types
        elif analysis_type == 'cdf-respt':
            if hasattr(solver, 'get_cdf_respt'):
                return solver.get_cdf_respt()
            return None

        elif analysis_type == 'cdf-passt':
            if hasattr(solver, 'get_cdf_passt'):
                return solver.get_cdf_passt()
            return None

        elif analysis_type == 'perct-respt':
            if hasattr(solver, 'get_perct_respt'):
                return solver.get_perct_respt(percentiles)
            return None

        # Transient analysis types
        elif analysis_type == 'tran-avg':
            if hasattr(solver, 'get_tran_avg'):
                solver.get_tran_avg()
                return "Transient analysis completed"
            return None

        elif analysis_type == 'tran-cdf-respt':
            if hasattr(solver, 'get_tran_cdf_respt'):
                return solver.get_tran_cdf_respt()
            return None

        elif analysis_type == 'tran-cdf-passt':
            if hasattr(solver, 'get_tran_cdf_passt'):
                return solver.get_tran_cdf_passt()
            return None

        # Probability analysis types
        elif analysis_type == 'prob':
            if hasattr(solver, 'get_prob') and node_idx is not None:
                if state:
                    return solver.get_prob(node_idx, state)
                else:
                    return solver.get_prob(node_idx)
            return None

        elif analysis_type == 'prob-aggr':
            if hasattr(solver, 'get_prob_aggr') and node_idx is not None:
                if state:
                    return solver.get_prob_aggr(node_idx, state)
                else:
                    return solver.get_prob_aggr(node_idx)
            return None

        elif analysis_type == 'prob-marg':
            if hasattr(solver, 'get_prob_marg') and node_idx is not None and class_idx is not None:
                if state:
                    return solver.get_prob_marg(node_idx, class_idx, state)
                else:
                    return solver.get_prob_marg(node_idx, class_idx)
            return None

        elif analysis_type == 'prob-sys':
            if hasattr(solver, 'get_prob_sys'):
                return solver.get_prob_sys()
            return None

        elif analysis_type == 'prob-sys-aggr':
            if hasattr(solver, 'get_prob_sys_aggr'):
                return solver.get_prob_sys_aggr()
            return None

        # Sampling analysis types (SSA only)
        elif analysis_type == 'sample':
            if hasattr(solver, 'sample') and node_idx is not None and model is not None:
                nodes = model.get_stateful_nodes()
                if node_idx < len(nodes):
                    return solver.sample(nodes[node_idx], num_events)
            return None

        elif analysis_type == 'sample-aggr':
            if hasattr(solver, 'sample_aggr') and node_idx is not None and model is not None:
                nodes = model.get_stateful_nodes()
                if node_idx < len(nodes):
                    return solver.sample_aggr(nodes[node_idx], num_events)
            return None

        elif analysis_type == 'sample-sys':
            if hasattr(solver, 'sample_sys'):
                return solver.sample_sys(num_events)
            return None

        elif analysis_type == 'sample-sys-aggr':
            if hasattr(solver, 'sample_sys_aggr'):
                return solver.sample_sys_aggr(num_events)
            return None

        # Reward analysis types (CTMC only)
        elif analysis_type == 'reward':
            if hasattr(solver, 'get_reward_result'):
                return solver.get_reward_result()
            return None

        elif analysis_type == 'reward-steady':
            if hasattr(solver, 'get_reward_steady_state'):
                return solver.get_reward_steady_state()
            return None

        elif analysis_type == 'reward-value':
            if hasattr(solver, 'get_reward_value_function') and reward_name:
                return solver.get_reward_value_function(reward_name)
            return None

    except Exception as e:
        return f"Error: {str(e)}"

    return None


def get_solver_results(solver, analysis_str, model=None, node_idx=None, class_idx=None,
                       state_str=None, num_events=1000, percentiles_str='50,90,95,99',
                       reward_name=None):
    """Get results from solver based on analysis types (comma-separated)."""
    results = {}

    # Parse analysis types
    analysis_types = [t.strip() for t in analysis_str.split(',')]
    state = parse_state(state_str) if state_str else None
    percentiles = parse_percentiles(percentiles_str)

    # Execute each analysis type
    for analysis_type in analysis_types:
        result = execute_single_analysis(
            solver, model, analysis_type,
            node_idx, class_idx, state, num_events, percentiles, reward_name
        )
        if result is not None:
            results[analysis_type] = result

    return results


def format_readable(results):
    """Format results as readable text tables."""
    output = []

    for key, value in results.items():
        if value is None:
            continue

        output.append(f"=== {key.upper()} ===")

        # Handle nested results (e.g., from 'all')
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                if sub_value is not None:
                    output.append(f"--- {sub_key} ---")
                    output.append(_format_single_result(sub_value))
                    output.append("")
        else:
            output.append(_format_single_result(value))
            output.append("")

    return "\n".join(output)


def _format_single_result(value):
    """Format a single result value."""
    if value is None:
        return "(no result)"

    # DataFrame (pandas)
    if hasattr(value, 'to_string'):
        return value.to_string(index=False)

    # String results
    if isinstance(value, str):
        return value

    # List/dict results
    if isinstance(value, (list, dict)):
        return json.dumps(value, indent=2, default=str)

    # Other objects - try to get string representation
    return str(value)


def format_json(results):
    """Format results as JSON."""
    json_results = {}

    for key, value in results.items():
        if value is not None:
            json_results[key] = _convert_to_json_serializable(value)

    return json.dumps(json_results, indent=2, default=str)


def _convert_to_json_serializable(value):
    """Convert a value to JSON-serializable format."""
    if value is None:
        return None

    # DataFrame (pandas)
    if hasattr(value, 'to_dict'):
        return value.to_dict(orient='records')

    # Nested dict (e.g., from 'all')
    if isinstance(value, dict):
        return {k: _convert_to_json_serializable(v) for k, v in value.items()}

    # List
    if isinstance(value, list):
        return [_convert_to_json_serializable(v) for v in value]

    # Numpy array
    if hasattr(value, 'tolist'):
        return value.tolist()

    # String, number, bool - return as-is
    if isinstance(value, (str, int, float, bool)):
        return value

    # Other objects - convert to string
    return str(value)


def format_csv(results):
    """Format results as CSV."""
    output = StringIO()

    for key, value in results.items():
        if value is None:
            continue

        output.write(f"# {key}\n")

        # Handle nested dicts (e.g., from 'all')
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                if sub_value is not None:
                    output.write(f"## {sub_key}\n")
                    _write_csv_value(output, sub_value)
                    output.write("\n")
        else:
            _write_csv_value(output, value)
            output.write("\n")

    return output.getvalue()


def _write_csv_value(output, value):
    """Write a value to CSV output."""
    # DataFrame (pandas)
    if hasattr(value, 'to_csv'):
        value.to_csv(output, index=False)
    # String
    elif isinstance(value, str):
        output.write(value)
        output.write("\n")
    # List/dict - write as JSON
    elif isinstance(value, (list, dict)):
        output.write(json.dumps(value, default=str))
        output.write("\n")
    else:
        output.write(str(value))
        output.write("\n")


def format_pickle(results):
    """Serialize results as pickle bytes."""
    return pickle.dumps(results)


def format_mat(model, output_filename):
    """Export model to MATLAB .mat format."""
    m2m = M2M()
    m2m.LINE2MAT(model, output_filename)


def format_results(results, output_format, model=None, output_file=None):
    """Format and output results based on specified format."""
    if output_format == 'readable':
        return format_readable(results)

    elif output_format == 'json':
        return format_json(results)

    elif output_format == 'csv':
        return format_csv(results)

    elif output_format == 'pickle':
        data = format_pickle(results)
        if output_file:
            with open(output_file, 'wb') as f:
                f.write(data)
            return f"Results saved to {output_file}"
        else:
            return data.hex()  # Return hex for stdout

    elif output_format == 'mat':
        if model is None:
            raise ValueError("Model required for .mat export")
        if output_file is None:
            raise ValueError("Output file required for .mat export")
        format_mat(model, output_file)
        return f"Model exported to {output_file}"

    else:
        raise ValueError(f"Unknown output format: {output_format}")


def run_cli(args=None):
    """Main CLI function."""
    parser = create_parser()
    args = parser.parse_args(args)

    # Set verbosity
    verbose = args.verbosity == 'normal'
    if verbose:
        GlobalConstants.setVerbose(VerboseLevel.STD)

    # Print banner if verbose
    if verbose:
        print("=" * 80)
        print("LINE Solver - Command Line Interface")
        print("=" * 80)
        print()

    try:
        # Server mode
        if args.port is not None:
            run_server(args.port, args.input)
            return

        # Validate analysis types
        analysis_types = validate_analysis_types(args.analysis)

        # Validate analysis-solver compatibility
        validate_analysis_solver_compat(analysis_types, args.solver)

        # Validate analysis parameters
        validate_analysis_params(analysis_types, args.node, args.class_idx, args.reward_name)

        # Standard CLI mode
        # Load model
        if args.file:
            model = load_model(args.file, args.input, verbose)
        else:
            # Read from stdin
            model = load_from_stdin(args.input, verbose)

        # Validate solver compatibility with input format
        input_format = args.input or detect_input_format(args.file)
        validate_solver_compatibility(input_format, args.solver)

        # Solve model
        solver = solve_model(model, args.solver, args.seed, verbose)

        # Get results with all new parameters
        results = get_solver_results(
            solver,
            args.analysis,
            model=model,
            node_idx=args.node,
            class_idx=args.class_idx,
            state_str=args.state,
            num_events=args.events,
            percentiles_str=args.percentiles,
            reward_name=args.reward_name
        )

        # Format and output results
        formatted = format_results(results, args.output, model)

        if isinstance(formatted, bytes):
            sys.stdout.buffer.write(formatted)
        else:
            print(formatted)

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def run_server(port=5863, default_input_format=None):
    """Run WebSocket server for remote model solving."""
    try:
        import websockets
    except ImportError:
        print("Error: websockets library required for server mode", file=sys.stderr)
        print("Install with: pip install websockets", file=sys.stderr)
        sys.exit(1)

    print("=" * 80)
    print(f"LINE Solver - Server Mode")
    print("=" * 80)
    print(f"Listening on port {port}")
    print("Press Ctrl+C to stop the server")
    print()

    async def handle_client(websocket, path):
        """Handle incoming WebSocket connection."""
        try:
            message = await websocket.recv()

            # Parse message: first line = args, rest = model content
            lines = message.split('\n', 1)
            args_str = lines[0]
            model_content = lines[1] if len(lines) > 1 else ""

            # Save model to temp file
            suffix = f'.{default_input_format or "jsim"}'
            with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as f:
                f.write(model_content)
                temp_filename = f.name

            try:
                # Parse arguments (comma-separated)
                arg_list = args_str.split(',')
                arg_list.extend(['--file', temp_filename])

                # Run CLI with these arguments
                parser = create_parser()
                parsed_args = parser.parse_args(arg_list)

                # Load and solve
                model = load_model(temp_filename, parsed_args.input, False)
                input_format = parsed_args.input or detect_input_format(temp_filename)
                validate_solver_compatibility(input_format, parsed_args.solver)

                solver = solve_model(model, parsed_args.solver, parsed_args.seed, False)
                results = get_solver_results(solver, parsed_args.analysis)
                formatted = format_results(results, parsed_args.output, model)

                # Send results back
                if isinstance(formatted, bytes):
                    await websocket.send(formatted.hex())
                else:
                    await websocket.send(formatted)

            finally:
                # Clean up temp file
                try:
                    os.unlink(temp_filename)
                except:
                    pass

        except Exception as e:
            print(f"Error handling client: {str(e)}", file=sys.stderr)
            try:
                await websocket.send(f"Error: {str(e)}")
            except:
                pass

    # Start server
    start_server = websockets.serve(handle_client, "0.0.0.0", port)

    loop = asyncio.get_event_loop()
    loop.run_until_complete(start_server)

    try:
        # Wait for Ctrl+C
        loop.run_forever()
    except KeyboardInterrupt:
        print("\nShutting down...")
    finally:
        loop.close()


class LineDockerClient:
    """Docker client for remote LINE solver communication.

    Sends models to a remote Docker-based LINE server for solving.
    Uses WebSocket protocol to communicate with the remote server.
    """

    def __init__(self, host='127.0.0.1', port='5863'):
        """Initialize Docker client with server connection details.

        Args:
            host (str): Docker server host (default: 127.0.0.1)
            port (str): Docker server port (default: 5863)
        """
        self.host = host
        self.port = port

    async def send_model(self, model_content, args=None, input_format='lqnx'):
        """Send model to remote Docker server and get results.

        Args:
            model_content (str): Model file content
            args (list): Command-line arguments for solver (e.g., ['-s', 'mva', '-o', 'json'])
            input_format (str): Input format (default: lqnx)

        Returns:
            str: Results from remote solver

        Raises:
            Exception: If connection or remote execution fails
        """
        try:
            import websockets
        except ImportError:
            raise ImportError("websockets library required. Install with: pip install websockets")

        # Build argument string
        if args is None:
            args = []
        args_str = ','.join(args) if args else ''

        # Build message: first line = args, rest = model content
        message = args_str + '\n' + model_content

        # Connect to remote server and send model
        uri = f"ws://{self.host}:{self.port}"

        try:
            async with websockets.connect(uri) as websocket:
                await websocket.send(message)
                result = await websocket.recv()
                return result
        except Exception as e:
            raise RuntimeError(f"Failed to connect to Docker server {uri}: {str(e)}")

    def send_model_sync(self, model_content, args=None, input_format='lqnx'):
        """Synchronous wrapper for send_model.

        Args:
            model_content (str): Model file content
            args (list): Command-line arguments
            input_format (str): Input format

        Returns:
            str: Results from remote solver
        """
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        try:
            return loop.run_until_complete(self.send_model(model_content, args, input_format))
        finally:
            loop.close()


def main():
    """Entry point for the CLI."""
    run_cli()


if __name__ == '__main__':
    main()
