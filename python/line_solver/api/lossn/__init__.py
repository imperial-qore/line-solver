"""
Loss Network Analysis Algorithms.

Native Python implementations for analyzing loss networks
using Erlang formulas and related methods.

Key algorithms:
    lossn_erlangfp: Erlang fixed-point algorithm for loss networks
    erlang_b: Erlang B blocking probability
    erlang_c: Erlang C delay probability
"""

from .erlang import (
    lossn_erlangfp,
    erlang_b,
    erlang_c,
)

__all__ = [
    'lossn_erlangfp',
    'erlang_b',
    'erlang_c',
]
