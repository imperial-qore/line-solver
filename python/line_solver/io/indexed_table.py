"""Re-export of the one IndexedTable, which lives at ``line_solver.indexed_table``.

This module used to hold a SECOND, DIVERGENT copy of the class: 329 lines
against the canonical 433, missing ``_matlab_format_float`` (so a result table
printed through it did not carry MATLAB's number formatting) and carrying its
own ``__getattr__``/``__iter__``. Nothing imported it -- every solver reaches
for ``line_solver.indexed_table`` -- but ``line_solver.io.__init__`` re-exports
the name, so ``line_solver.io.IndexedTable`` was a public class with a second
implementation behind it. Whichever of the two a caller happened to get decided
how its numbers printed.

The MATLAB tree carried the same pair (``lang/IndexedTable.m`` beside
``io/IndexedTable.m``), byte-identical there and removed on 2026-08-16; this one
had drifted, so it is collapsed rather than deleted, and the public name keeps
working.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from ..indexed_table import IndexedTable

__all__ = ['IndexedTable']
