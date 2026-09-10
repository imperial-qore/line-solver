"""
There must be exactly ONE TimingStrategy enum in native python.

Two classes of the same name used to exist, with DIFFERENT values:
line_solver.lang.nodes as an IntEnum (TIMED = 0, IMMEDIATE = 1, matching MATLAB
and the JAR) and line_solver.constants through auto() (TIMED = 1, IMMEDIATE = 2).
The package root exported the first, while io/__init__.py's M2M.JSIM2LINE
imported the second, so a JSIM-imported immediate transition carried a value that
never compared equal to TimingStrategy.IMMEDIATE anywhere downstream and was
served as a TIMED transition. Neither class is comparable to the other and both
print the same repr, so the mismatch is invisible in a debugger.

This guards the shape of the defect rather than one of its symptoms: a second
declaration reintroduced anywhere in the package fails here.
"""

import os
import pkgutil
import importlib
import inspect

import pytest

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import TimingStrategy
from line_solver.lang.nodes import TimingStrategy as NodesTimingStrategy


def test_package_root_exports_the_model_layer_enum():
    assert TimingStrategy is NodesTimingStrategy


def test_values_match_matlab_and_the_jar():
    assert int(TimingStrategy.TIMED) == 0
    assert int(TimingStrategy.IMMEDIATE) == 1


def test_constants_declares_no_second_timing_strategy():
    import line_solver.constants as constants
    assert not hasattr(constants, 'TimingStrategy'), (
        'constants.py declares a second TimingStrategy again; the model layer '
        'compares against line_solver.lang.nodes.TimingStrategy')


def test_exactly_one_timing_strategy_class_in_the_package():
    found = {}
    root = os.path.dirname(line_solver.__file__)
    for mod in pkgutil.walk_packages([root], prefix='line_solver.'):
        name = mod.name
        # Third-party vendored trees and optional backends are not ours to police.
        if '.lib.' in name or name.endswith('.lib'):
            continue
        try:
            module = importlib.import_module(name)
        except Exception:
            continue
        obj = getattr(module, 'TimingStrategy', None)
        if obj is None or not inspect.isclass(obj):
            continue
        # Only the module that DEFINES it counts; every other module re-exports.
        if getattr(obj, '__module__', None) == name:
            found[name] = obj
    assert list(found) == ['line_solver.lang.nodes'], (
        'TimingStrategy is defined in %s; exactly one definition is allowed' % sorted(found))
