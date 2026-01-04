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
        if options is not None and 'verbose' in options:
            opt_verbose = options['verbose']
            if isinstance(opt_verbose, VerboseLevel):
                is_debug = is_debug or (opt_verbose == VerboseLevel.DEBUG)
            elif isinstance(opt_verbose, int):
                is_debug = is_debug or (opt_verbose == VerboseLevel.DEBUG.value)

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
    'line_verbosity',
    'get_logger',
]
