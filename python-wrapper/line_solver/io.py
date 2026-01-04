"""
LINE Solver I/O Functions

This module provides I/O functions for exporting LINE network models to various formats.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import jpype


def QN2JSIMG(model, outputFileName=None):
    """
    Writes a Network model to JMT JSIMG format.

    Args:
        model: Network model to export
        outputFileName: Optional output file path (default: temp file)

    Returns:
        Path to the generated JSIMG file

    Example:
        >>> model = Network('example')
        >>> # ... define model ...
        >>> fname = QN2JSIMG(model)
        >>> print(f"Model written to: {fname}")
    """
    # Get the Java io.QN2JSIMG class
    java_qn2jsimg = jpype.JPackage('jline').io.QN2JSIMG

    # Get the underlying Java model object
    java_model = model.obj if hasattr(model, 'obj') else model

    if outputFileName is None:
        return java_qn2jsimg.writeJSIM(java_model)
    else:
        return java_qn2jsimg.writeJSIM(java_model, outputFileName)


# Alias for snake_case
qn2jsimg = QN2JSIMG


class JLINE:
    """
    Utility class for converting between Python LINE models and Java JLINE representations.

    This class provides static methods for bidirectional conversion, enabling
    interoperability with Java-based solvers and tools.
    """

    @staticmethod
    def line_to_jline(model):
        """
        Convert a Python Network model to its Java JLINE representation.

        Args:
            model (Network): Python Network model to convert.

        Returns:
            jline.lang.Network: The underlying Java Network object.

        Example:
            >>> model = Network('Example')
            >>> java_model = JLINE.line_to_jline(model)
        """
        return model.obj

    @staticmethod
    def jline_to_line(jnetwork):
        """
        Convert a Java JLINE Network to a Python Network wrapper.

        Args:
            jnetwork: Java jline.lang.Network object.

        Returns:
            Network: Python Network wrapper around the Java object.

        Example:
            >>> recovered_model = JLINE.jline_to_line(java_model)
        """
        from .lang import Network
        # Handle wrapped JNetwork if needed
        if hasattr(jnetwork, 'obj'):
            jnetwork = jnetwork.obj
        return Network(jnetwork)


def LQN2QN(lqn):
    """
    Convert a LayeredNetwork (LQN) to a Network (QN).

    This converter transforms a Layered Queueing Network into an equivalent
    Queueing Network that approximates LQN blocking semantics using a blocking
    approximation technique.

    Blocking Approximation:
    The conversion models synchronous call blocking by aggregating service
    times: when a task makes a synchCall, its effective service time includes
    both its own processing and all downstream processing. This approach:
    - Correctly models throughput (validated against LQNS)
    - Approximates response times (aggregated into entry task)
    - Avoids explicit queueing constraints

    Phase-2 Support:
    When an entry has phase-2 activities, a sequential Fork-Join structure
    is created to model the asynchronous phase-2 processing, allowing work
    to continue after sending a reply back to the caller.

    Supported Features:
    - Multi-tier architectures (client-server, 3-tier, etc.)
    - Synchronous calls (synchCall) with nested depth
    - Phase-2 activities (asynchronous work after reply)
    - Multiple processors and multiplicity
    - Different scheduling disciplines (FCFS, etc.)

    Limitations:
    - Approximation becomes less accurate for complex blocking patterns
    - Does not model priorities or preemption
    - Response time aggregation may be conservative

    Args:
        lqn: LayeredNetwork model to convert.

    Returns:
        Network: Equivalent queueing network approximating the LQN behavior.

    Raises:
        Exception: If the LQN structure is invalid or conversion fails.

    Examples:
        >>> from line_solver import *
        >>> lqn = LayeredNetwork('MyLQN')
        >>> # ... define 2-tier LQN model ...
        >>> model = lqn2qn(lqn)  # Use snake_case alias
        >>> solver = SolverDES(model)
        >>> results = solver.getAvgTable()

        >>> # For 3-tier application
        >>> lqn = LayeredNetwork('AppLQN')
        >>> # ... define 3-tier model (web, app, database) ...
        >>> model = LQN2QN(lqn)
        >>> solver = SolverMVA(model)
        >>> results = solver.getAvgTable()

    Note:
        The converter validates LQN synchronous call chains and automatically
        aggregates service times. For complex models, use SolverDES or SolverJMT
        with the converted model for better accuracy.

        See also: example_lqn2qn_basic.py for complete working examples.
    """
    from .lang import Network

    # Get the underlying Java LQN object
    java_lqn = lqn.obj if hasattr(lqn, 'obj') else lqn

    # Call Java LQN2QN converter
    java_qn = jpype.JPackage('jline').io.LQN2QN.convert(java_lqn)

    # Wrap the result in a Python Network
    return Network(java_qn)


# Alias for snake_case
lqn2qn = LQN2QN
