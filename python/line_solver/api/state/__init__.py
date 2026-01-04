"""
State analysis framework for LINE networks (pure Python).

This package provides state space analysis and probability computation
for queueing networks.
"""

from .marginal import toMarginal, fromMarginal
from .space_generator import spaceGenerator

__all__ = [
    'toMarginal',
    'fromMarginal',
    'spaceGenerator',
]
