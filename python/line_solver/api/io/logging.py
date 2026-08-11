"""
LINE logging and output utilities.

This module provides consistent logging and output functions for the LINE solver,
matching the MATLAB io/ utilities: line_warning, line_error, line_printf,
line_debug, and line_verbosity.

Port from:
    - matlab/src/io/line_warning.m
    - matlab/src/io/line_error.m
    - matlab/src/io/line_printf.m
    - matlab/src/io/line_debug.m
    - matlab/src/io/line_verbosity.m
"""

import sys
import time
import warnings
import traceback
from enum import Enum
from typing import Optional, Any, Dict
from dataclasses import dataclass, field


class VerboseLevel(Enum):
    """Verbosity levels for LINE output."""
    SILENT = 0
    STD = 1
    DEBUG = 2


@dataclass
class _WarningState:
    """State for warning suppression."""
    last_warning: str = ''
    suppressed_warnings: bool = False
    suppressed_announcement: bool = False
    last_warning_time: float = 0.0
    suppression_start_time: float = 0.0


class LineLogger:
    """
    Singleton logger for LINE solver output.

    Provides consistent logging, warning, and error handling across the
    LINE Python native implementation.

    Attributes:
        verbose: Current verbosity level
        stdout: Output stream (default sys.stdout)
        warning_state: State for warning suppression
    """

    _instance: Optional['LineLogger'] = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if self._initialized:
            return
        self._initialized = True
        self.verbose = VerboseLevel.STD
        self.stdout = sys.stdout
        self._warning_state = _WarningState()
        self._suppression_timeout = 60.0  # seconds

    @classmethod
    def get_instance(cls) -> 'LineLogger':
        """Get the singleton logger instance."""
        return cls()

    def set_verbosity(self, level: VerboseLevel) -> None:
        """
        Set the verbosity level.

        Args:
            level: VerboseLevel enum value
        """
        self.verbose = level

        # Configure Python warnings based on verbosity
        if level == VerboseLevel.SILENT:
            warnings.filterwarnings('ignore')
        else:
            warnings.filterwarnings('default')

    def get_verbosity(self) -> VerboseLevel:
        """Get the current verbosity level."""
        return self.verbose

    def printf(self, msg: str, *args, **kwargs) -> None:
        """
        Print formatted message if not in SILENT mode.

        Args:
            msg: Format string
            *args: Format arguments
            **kwargs: Additional keyword arguments (ignored)
        """
        if self.verbose == VerboseLevel.SILENT:
            return

        if args:
            try:
                formatted = msg % args
            except (TypeError, ValueError):
                formatted = msg.format(*args) if '{' in msg else msg
        else:
            formatted = msg

        self.stdout.write(formatted)
        self.stdout.flush()

    def warning(self, caller: str, msg: str, *args) -> None:
        """
        Print warning message with caller information.

        Implements warning suppression to avoid flooding output with
        repeated warnings. Same warning is suppressed for 60 seconds.

        Args:
            caller: Name of the calling function/module
            msg: Warning message format string
            *args: Format arguments
        """
        if self.verbose == VerboseLevel.SILENT:
            return

        # Format the message
        if args:
            try:
                formatted_msg = msg % args
            except (TypeError, ValueError):
                formatted_msg = msg.format(*args) if '{' in msg else msg
        else:
            formatted_msg = msg

        final_msg = f"Warning [{caller}]: {formatted_msg}"
        current_time = time.time()

        # Check if this is a repeated warning
        state = self._warning_state
        time_since_suppression = current_time - state.suppression_start_time
        time_since_last = current_time - state.last_warning_time

        if final_msg != state.last_warning or time_since_suppression > self._suppression_timeout:
            # New warning or timeout expired
            self.printf(f"{final_msg}\n")
            state.last_warning = final_msg
            state.suppressed_warnings = False
            state.suppressed_announcement = False
            state.suppression_start_time = current_time
        else:
            # Repeated warning
            if not state.suppressed_warnings and not state.suppressed_announcement:
                self.printf(
                    f"[{caller}] Message cast more than once, "
                    f"repetitions will not be printed for {int(self._suppression_timeout)} seconds.\n"
                )
                state.suppressed_announcement = True
                state.suppressed_warnings = True
                state.suppression_start_time = current_time

        state.last_warning_time = current_time

    def warning_always(self, caller: str, msg: str, *args) -> None:
        """
        Print a warning without the repeat suppression applied by warning().

        warning() keeps only the last message and hides an identical repeat for
        60 seconds. That is right for configuration notices cast once per model,
        but wrong for a warning that reports a correctness limitation of the
        analysis: solving several models in one session would then flag only the
        first one, and the user would read the silence on the others as a clean
        bill of health. Warnings that say "these numbers are not exact" must be
        raised for every model they apply to, so they go through here instead.
        Verbosity gating is unchanged: SILENT still silences everything.

        Args:
            caller: Name of the calling function/module
            msg: Warning message format string
            *args: Format arguments
        """
        if self.verbose == VerboseLevel.SILENT:
            return

        if args:
            try:
                formatted_msg = msg % args
            except (TypeError, ValueError):
                formatted_msg = msg.format(*args) if '{' in msg else msg
        else:
            formatted_msg = msg

        self.printf(f"Warning [{caller}]: {formatted_msg}\n")

    def error(self, caller: str, msg: str) -> None:
        """
        Raise an error with caller information.

        Args:
            caller: Name of the calling function/module
            msg: Error message

        Raises:
            RuntimeError: Always raises with formatted message
        """
        # Get call stack info
        stack = traceback.extract_stack()
        if len(stack) >= 2:
            line_num = stack[-2].lineno
            filename = stack[-2].filename
        else:
            line_num = 0
            filename = caller

        error_str = f"[{caller} @ line {line_num}] {msg}"
        raise RuntimeError(error_str)

    def debug(self, msg: str, *args, options: Optional[Dict] = None) -> None:
        """
        Print debug message if in DEBUG mode.

        Args:
            msg: Debug message format string
            *args: Format arguments
            options: Optional dict with 'verbose' key to override global setting
        """
        # Check if debug mode is enabled
        is_debug = self.verbose == VerboseLevel.DEBUG
        if options is not None:
            # Support both dict-like options and dataclass/object options
            opt_verbose = None
            if isinstance(options, dict):
                opt_verbose = options.get('verbose')
            elif hasattr(options, 'verbose'):
                opt_verbose = options.verbose
            if opt_verbose is not None:
                if isinstance(opt_verbose, VerboseLevel):
                    is_debug = is_debug or (opt_verbose == VerboseLevel.DEBUG)
                elif isinstance(opt_verbose, (int, bool)):
                    is_debug = is_debug or (int(opt_verbose) == VerboseLevel.DEBUG.value)

        if not is_debug:
            return

        # Format the message
        if args:
            try:
                formatted_msg = msg % args
            except (TypeError, ValueError):
                formatted_msg = msg.format(*args) if '{' in msg else msg
        else:
            formatted_msg = msg

        self.printf(f"[DEBUG] {formatted_msg}\n")


# Global logger instance
_logger = LineLogger()


# Module-level convenience functions

def line_printf(msg: str, *args) -> None:
    """
    Print formatted message if not in SILENT mode.

    Args:
        msg: Format string
        *args: Format arguments

    References:
        MATLAB: matlab/src/io/line_printf.m
    """
    _logger.printf(msg, *args)


def line_warning(caller: str, msg: str, *args) -> None:
    """
    Print warning message with caller information.

    Implements warning suppression to avoid repeated warnings.

    Args:
        caller: Name of the calling function/module
        msg: Warning message format string
        *args: Format arguments

    References:
        MATLAB: matlab/src/io/line_warning.m
    """
    _logger.warning(caller, msg, *args)


def line_warning_always(caller: str, msg: str, *args) -> None:
    """
    Print a warning without line_warning's repeat suppression.

    Use for warnings that report a correctness limitation of the analysis, which
    must be raised for every model they apply to rather than be hidden as a
    repeat of the previous model's warning.

    Args:
        caller: Name of the calling function/module
        msg: Warning message format string
        *args: Format arguments

    References:
        MATLAB: matlab/src/io/line_warning_always.m
    """
    _logger.warning_always(caller, msg, *args)


def line_error(caller: str, msg: str) -> None:
    """
    Raise an error with caller information.

    Args:
        caller: Name of the calling function/module
        msg: Error message

    Raises:
        RuntimeError: Always raises with formatted message

    References:
        MATLAB: matlab/src/io/line_error.m
    """
    _logger.error(caller, msg)


def line_debug(msg: str, *args, options: Optional[Dict] = None) -> None:
    """
    Print debug message if in DEBUG mode.

    Args:
        msg: Debug message format string
        *args: Format arguments
        options: Optional dict with 'verbose' key

    References:
        MATLAB: matlab/src/io/line_debug.m
    """
    _logger.debug(msg, *args, options=options)


# Attribution strings for the external tools the in-tree wrapper solvers
# delegate to. Verified against each upstream project's own pages; do not reword
# an author list or a URL from memory. Mirror any edit in the MATLAB
# (matlab/src/io/line_ack.m) and Java (jline.io.InputOutput.line_ack) tables.
# Out-of-tree wrapper solvers pass their own text through the msg argument:
# their tools must not be named in this codebase.
_TOOL_ACK = {
    'JMT': (
        "SolverJMT delegates to Java Modelling Tools (JMT), by M. Bertoli, "
        "G. Casale, G. Serazzi (Politecnico di Milano, Imperial College London). "
        "Please acknowledge the JMT authors: http://jmt.sourceforge.net/"
    ),
    'LQNS': (
        "SolverLQNS delegates to LQNS/LQSIM, by G. Franks, M. Woodside et al. "
        "(Real-Time and Distributed Systems Group, Carleton University). "
        "Please acknowledge the LQNS authors: http://www.layeredqueues.org/"
    ),
    'QNS': (
        "SolverQNS delegates to qnsolver, part of the LQNS distribution by "
        "G. Franks, M. Woodside et al. (Real-Time and Distributed Systems Group, "
        "Carleton University). Please acknowledge the LQNS authors: "
        "http://www.layeredqueues.org/"
    ),
}

# One-line reference to the canonical paper of each tool, printed under the
# acknowledgement. It describes the same work as the BibTeX entry in
# _TOOL_BIBTEX (keys BerCS07 and fran.ea09 in doc/latex/biblio.bib), which is
# where the citation key belongs: the printed line is for the reader.
_TOOL_CITE = {
    'JMT': (
        'M. Bertoli, G. Casale, G. Serazzi. "The JMT Simulator for Performance '
        'Evaluation of Non-Product-Form Queueing Networks". Proc. of the 40th '
        'Annual Simulation Symposium (ANSS), pp. 3-10, 2007.'
    ),
    'LQNS': (
        'G. Franks, T. Al-Omari, M. Woodside, O. Das, S. Derisavi. "Enhanced '
        'Modeling and Solution of Layered Queueing Networks". IEEE Trans. '
        'Software Engineering, 35(2):148-161, 2009.'
    ),
}
# qnsolver ships inside the LQNS distribution and shares its reference.
_TOOL_CITE['QNS'] = _TOOL_CITE['LQNS']

# Machine-readable form of the same references, so an acknowledgement can be
# turned into a citation without retyping it. Copied verbatim from
# doc/latex/biblio.bib; mirror any edit in matlab/src/io/line_citation.m and
# jline.io.InputOutput.line_citation.
_TOOL_BIBTEX = {
    'JMT': (
        "@INPROCEEDINGS{BerCS07,\n"
        "  author = {M. Bertoli and G. Casale and G. Serazzi},\n"
        "  title = {The {JMT} Simulator for Performance Evaluation of Non-Product-Form\n"
        "\tQueueing Networks},\n"
        "  booktitle = {Proc. of the 40th Annual Simulation Symposium (ANSS)},\n"
        "  year = {2007},\n"
        "  pages = {3--10}\n"
        "}"
    ),
    'LQNS': (
        "@ARTICLE{fran.ea09,\n"
        "  author = {G. Franks and T. Al-Omari and M. Woodside and O. Das and S. Derisavi},\n"
        "  title = {Enhanced Modeling and Solution of Layered Queueing Networks},\n"
        "  journal = {IEEE Trans. Software Engineering},\n"
        "  year = {2009},\n"
        "  volume = {35},\n"
        "  pages = {148-161},\n"
        "  number = {2}\n"
        "}"
    ),
}
_TOOL_BIBTEX['QNS'] = _TOOL_BIBTEX['LQNS']

# Tools already acknowledged in this session (one interpreter = one session).
# Distinct from a collective flag: each external tool is acknowledged on its own.
_tool_ack_shown = set()


def _ack_is_debug(verbose: Any) -> bool:
    """Test whether the caller asked for DEBUG verbosity, across its encodings.

    Wrapper options carry verbosity as a bool (most solvers), as the string
    'silent'/'std'/'debug' (LDES), or as a VerboseLevel; None means "inherit
    the global". A bool can only express SILENT or STD, so it never selects
    DEBUG. The VerboseLevel test goes by name because line_solver.constants
    defines a second enum of the same name.
    """
    if verbose is None:
        return _logger.get_verbosity() == VerboseLevel.DEBUG
    if isinstance(verbose, bool):
        return False
    if isinstance(verbose, str):
        return verbose.lower() == 'debug'
    if isinstance(verbose, Enum):
        return verbose.name == 'DEBUG'
    return verbose == VerboseLevel.DEBUG.value


def line_ack(tool_name: str, verbose: Any = None, msg: Optional[str] = None,
             cite: Optional[str] = None) -> None:
    """
    Print, once per session, the acknowledgement of the external tool that a
    wrapper solver delegates to, together with the pointer to its official
    website and the canonical paper to cite.

    The acknowledgement is pull-based, like the library attribution: nothing is
    printed at the default verbosity, and only a caller that asks for it
    explicitly, by running at VerboseLevel.DEBUG, gets the line, at most once
    per tool per session. solver.citations() and line_citation() are the quiet
    ways to obtain the same reference.

    Args:
        tool_name: External tool, e.g. 'JMT', 'LQNS', 'QNS'
        verbose: Caller's verbosity (bool, str or VerboseLevel); None inherits
            the global level
        msg: Acknowledgement text for a solver living outside this tree; when
            omitted the text comes from the in-tree table
        cite: One-line reference for a solver living outside this tree; when
            omitted it comes from the in-tree table

    References:
        MATLAB: matlab/src/io/line_ack.m
    """
    if not _ack_is_debug(verbose) or _logger.get_verbosity() == VerboseLevel.SILENT:
        return
    key = tool_name.upper()
    if msg is None:
        msg = _TOOL_ACK.get(key)
    if cite is None:
        cite = _TOOL_CITE.get(key)
    if msg is None or key in _tool_ack_shown:
        return
    _tool_ack_shown.add(key)
    _logger.printf('%s\n', msg)
    if cite:
        _logger.printf('  Cite: %s\n', cite)


def line_citation(tool_name: str) -> str:
    """
    Return the BibTeX entry for the canonical paper of an external tool that a
    wrapper solver delegates to, so that the acknowledgement printed by
    line_ack can be turned into a citation without retyping it.

    Args:
        tool_name: External tool, 'JMT', 'LQNS' or 'QNS' ('QNS' shares the LQNS
            reference, qnsolver being part of that distribution)

    Returns:
        The BibTeX entry, or '' for an unknown tool. Wrapper solvers that live
        outside this tree carry their own reference: their tools must not be
        named in this codebase.

    References:
        MATLAB: matlab/src/io/line_citation.m
    """
    return _TOOL_BIBTEX.get(tool_name.upper(), '')


def line_verbosity(level: VerboseLevel = VerboseLevel.STD) -> None:
    """
    Set the global verbosity level for LINE.

    Args:
        level: VerboseLevel enum value (default: STD)

    References:
        MATLAB: matlab/src/io/line_verbosity.m
    """
    _logger.set_verbosity(level)


def get_logger() -> LineLogger:
    """Get the global LINE logger instance."""
    return _logger


__all__ = [
    'VerboseLevel',
    'LineLogger',
    'line_printf',
    'line_warning',
    'line_error',
    'line_debug',
    'line_ack',
    'line_citation',
    'line_verbosity',
    'get_logger',
]
